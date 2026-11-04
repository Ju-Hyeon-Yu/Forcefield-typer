// 웹 타이퍼를 Python 결과와 대조하고 입력 형식·3D 생성 품질을 확인하는 Node 테스트
const fs = require("fs");
const os = require("os");
const path = require("path");
const assert = require("assert");
const { execFileSync } = require("child_process");
const T = require("../src/typer.js");
const IO = require("../src/io.js");
const P = require("../src/pipeline.js");
const G2 = require("../src/gaff2.js");
const initRDKit = require("@rdkit/rdkit");

const HERE = __dirname;
const TABLES = path.join(HERE, "../build/tables.json");
const GAFF2_TABLES = path.join(HERE, "../build/gaff2_tables.json");
const GAFF2_SMILES = path.join(HERE, "../../tests/data/gaff2_smiles.txt");
const results = [];
const test = async (name, fn) => {
  try { await fn(); results.push([true, name]); console.log("ok   ", name); }
  catch (e) { results.push([false, name]); console.log("FAIL ", name, "\n     ", e.message.split("\n").slice(0, 8).join("\n      ")); }
};

(async () => {
  if (!fs.existsSync(TABLES)) {
    fs.mkdirSync(path.dirname(TABLES), { recursive: true });
    fs.writeFileSync(TABLES, execFileSync("python3", [path.join(HERE, "../build_tables.py")]));
  }
  if (!fs.existsSync(GAFF2_TABLES)) {
    fs.writeFileSync(GAFF2_TABLES, execFileSync("python3", [path.join(HERE, "../build_tables.py"), "--ff", "gaff2"]));
  }
  const R = await initRDKit();
  const ff = new T.DreidingFF(JSON.parse(fs.readFileSync(TABLES)));
  const ffs = { dreiding: ff, gaff2: new G2.Gaff2(JSON.parse(fs.readFileSync(GAFF2_TABLES))) };
  const tmp = fs.mkdtempSync(path.join(os.tmpdir(), "dreiding-web-"));

  await test("라벨·경고·오류가 Python 타이퍼와 같음 (SMILES 142개)", () => {
    const py = JSON.parse(execFileSync("python3", [path.join(HERE, "py_labels.py"), path.join(HERE, "smiles_set.txt")]));
    const js = JSON.parse(execFileSync("node", [path.join(HERE, "js_labels.js"), path.join(HERE, "smiles_set.txt")]));
    const diff = Object.keys(py).filter(k => JSON.stringify(py[k]) !== JSON.stringify(js[k]));
    assert.strictEqual(diff.length, 0, "다른 SMILES: " + diff.join(" "));
  });

  await test("GAFF2 타입·경고·보충 블록·오류가 Python 과 같음 (antechamber 정답 분자 391개)", () => {
    const py = JSON.parse(execFileSync("python3", [path.join(HERE, "py_gaff2.py"), GAFF2_SMILES],
      { maxBuffer: 1 << 26, stdio: ["ignore", "pipe", "ignore"] }));
    const js = JSON.parse(execFileSync("node", [path.join(HERE, "js_gaff2.js"), GAFF2_SMILES], { maxBuffer: 1 << 26 }));
    assert.ok(Object.keys(py).length > 380);
    const diff = Object.keys(py).filter(k => JSON.stringify(py[k]) !== JSON.stringify(js[k]));
    assert.strictEqual(diff.length, 0, "다른 SMILES: " + diff.join(" "));
  });

  execFileSync("python3", [path.join(HERE, "make_cases.py"), tmp]);
  const cases = JSON.parse(fs.readFileSync(path.join(tmp, "cases.json")));
  const read = p => ({ name: path.basename(p), text: fs.readFileSync(p, "utf8") });
  const body = s => s.split("\n").slice(1).join("\n").replace(/^# 전하: .*$/m, "");

  await test(`mol2 입력 .lt 가 Python CLI 와 같음 (원자 순서 섞음, ${cases.length}개 계)`, () => {
    for (const c of cases) {
      const r = P.run(R, ff, { smiles: c.smiles, name: "MOL", file: read(c.mol2) });
      assert.strictEqual(body(r.lt), body(fs.readFileSync(c.lt, "utf8")), c.smiles);
    }
  });

  await test("GAFF2: mol2 입력 .lt(보충 블록 포함)가 Python CLI 와 같음", () => {
    let n = 0;
    for (const c of cases) {
      if (!c.lt_gaff2) {
        assert.throws(() => P.run(R, ffs, { smiles: c.smiles, name: "MOL", file: read(c.mol2), ff: "gaff2" }), T.TypingError);
        continue;
      }
      const r = P.run(R, ffs, { smiles: c.smiles, name: "MOL", file: read(c.mol2), ff: "gaff2" });
      assert.strictEqual(body(r.lt), body(fs.readFileSync(c.lt_gaff2, "utf8")), c.smiles);
      n++;
    }
    assert.ok(n >= cases.length - 2 && cases.some(c => /\nGAFF2 \{\n/.test(fs.readFileSync(c.lt_gaff2 || c.lt, "utf8"))));
  });

  await test("pdb(CONECT 있음/없음)·xyz 입력이 mol2 와 같은 라벨·좌표 (전하 0)", () => {
    for (const c of cases) {
      const ref = P.run(R, ff, { smiles: c.smiles, file: read(c.mol2) });
      for (const key of ["pdb", "pdb_noconect", "xyz"]) {
        const r = P.run(R, ff, { smiles: c.smiles, file: read(c[key]) });
        assert.deepStrictEqual(r.labels, ref.labels, `${key} ${c.smiles}`);
        r.rows.forEach((row, i) => {
          assert.ok(row.xyz.every((v, k) => Math.abs(v - ref.rows[i].xyz[k]) < 1e-3), `${key} 좌표 ${c.smiles}`);
          assert.strictEqual(row.q, 0);
        });
        assert.ok(r.warnings.some(w => w.includes("전하 정보가 없어")));
      }
    }
  });

  await test("SMILES 와 좌표 파일 구조가 다르면 오류", () => {
    const c = cases[0];
    assert.throws(() => P.run(R, ff, { smiles: "CCO", file: read(c.mol2) }), IO.InputError);
    assert.throws(() => P.run(R, ff, { smiles: "c1ccccc1.[I-]", file: read(cases[7].mol2) }), IO.InputError);
  });

  await test("잘못된 입력은 알아들을 수 있는 오류", () => {
    assert.throws(() => P.run(R, ff, { smiles: "" }), /SMILES 를 입력/);
    assert.throws(() => P.run(R, ff, { smiles: "C1CC" }), /읽지 못했/);
    assert.throws(() => P.run(R, ff, { smiles: "CC", name: "1abc" }), /분자 이름/);
    assert.throws(() => P.run(R, ff, { smiles: "CC", file: { name: "x.sdf", text: "" } }), /지원하지 않는/);
    assert.throws(() => P.run(R, ff, { smiles: "F[P-](F)(F)(F)(F)F" }), T.TypingError);
  });

  await test("SMILES 만일 때 3D 구조 품질 (결합 길이·각·원자 겹침·조각 간격)", () => {
    const smis = fs.readFileSync(path.join(HERE, "smiles_set.txt"), "utf8").split(/\s+/).filter(Boolean);
    let worstBond = 0, worstAng = 0, minNb = 99, n = 0, t0 = Date.now(), worst = "";
    for (const smi of smis) {
      let r;
      try { r = P.run(R, ff, { smiles: smi }); } catch (e) { continue; }
      n++;
      const X = r.rows.map(x => x.xyz), g = r.graph;
      assert.ok(X.every(p => p.every(Number.isFinite)), "NaN " + smi);
      const d = (i, j) => Math.hypot(X[i][0] - X[j][0], X[i][1] - X[j][1], X[i][2] - X[j][2]);
      for (const b of g.bonds) {
        const dev = Math.abs(d(b.a, b.b) - ff.bondR0(r.labels[b.a], r.labels[b.b]));
        if (dev > worstBond) { worstBond = dev; worst = smi; }
      }
      const nbset = g.nbrs.map(x => new Set(x.map(([j]) => j)));
      const smallRing = c => [...nbset[c]].some(p => [...nbset[c]].some(q => p !== q &&
        (nbset[p].has(q) || [...nbset[p]].some(m => m !== c && nbset[q].has(m)))));
      g.atoms.forEach((a, c) => {   // 사면체 C 의 결합각 (3·4원 고리 제외)
        if (r.labels[c] !== "C_3" || g.deg[c] !== 4 || smallRing(c)) return;
        const nb = g.nbrs[c].map(([j]) => j);
        for (let p = 0; p < 4; p++) for (let q = p + 1; q < 4; q++) {
          const u = X[nb[p]].map((v, k) => v - X[c][k]), v = X[nb[q]].map((v2, k) => v2 - X[c][k]);
          const ang = Math.acos((u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (Math.hypot(...u) * Math.hypot(...v))) * 180 / Math.PI;
          worstAng = Math.max(worstAng, Math.abs(ang - 109.47));
        }
      });
      const topo = g.atoms.map((_, s) => { const dd = new Array(g.atoms.length).fill(99); dd[s] = 0; const q = [s];
        while (q.length) { const u = q.shift(); for (const [v] of g.nbrs[u]) if (dd[v] === 99) { dd[v] = dd[u] + 1; q.push(v); } } return dd; });
      for (let i = 0; i < X.length; i++) for (let j = i + 1; j < X.length; j++) {
        if (topo[i][j] >= 3) minNb = Math.min(minNb, d(i, j));
      }
    }
    console.log(`      분자 ${n}개, 최대 결합 편차 ${worstBond.toFixed(3)} Å (${worst}), 사면체 C 최대 각 편차 ` +
      `${worstAng.toFixed(1)}°, 최소 비결합(1-4 이상) 거리 ${minNb.toFixed(2)} Å, ${(Date.now() - t0)} ms`);
    assert.ok(worstBond < 0.05, "결합 길이");
    assert.ok(worstAng < 10, "사면체 각");
    assert.ok(minNb > 1.7, "원자 겹침");
  });

  await test("GAFF2: SMILES 만일 때도 3D 좌표를 만듦 (DREIDING 라벨이 안 되는 CO2·PF6⁻·아이오도늄 포함)", () => {
    const smis = ["O=C=O", "C=C=C", "F[P-](F)(F)(F)(F)F", "c1ccc(cc1)[I+]c1ccccc1", "[O-]Cl(=O)(=O)=O",
      "CN=C=O", "CC(=O)C=[N+]=[N-]", "c1ccc-2c(c1)-c1ccccc-21"];
    for (const smi of smis) {
      assert.throws(() => P.run(R, ffs, { smiles: smi }), T.TypingError, "DREIDING 은 실패해야 함 " + smi);
    }
    for (const smi of smis.concat(["CC(N)=O", "c1ccccc1"])) {
      const r = P.run(R, ffs, { smiles: smi, ff: "gaff2" });
      const X = r.rows.map(x => x.xyz);
      assert.ok(X.every(p => p.every(Number.isFinite)), "NaN " + smi);
      for (let i = 0; i < X.length; i++) for (let j = i + 1; j < X.length; j++) {
        const d = Math.hypot(X[i][0] - X[j][0], X[i][1] - X[j][1], X[i][2] - X[j][2]);
        assert.ok(d > 0.85, `원자 겹침 ${smi} ${i + 1}-${j + 1} ${d.toFixed(2)}`);
      }
      assert.ok(r.lt.includes('import "gaff2.lt"') && r.lt.includes("MOL inherits GAFF2"));
    }
  });

  await test("SMILES 만일 때 이온 조각 사이가 3 Å 이상", () => {
    const r = P.run(R, ff, { smiles: "Cn1cc[n+](Cc2ccc(I)cc2)c1.[O-]S(=O)(=O)C(F)(F)F.[I-]" });
    const fr = T.fragments(r.graph), X = r.rows.map(x => x.xyz);
    for (let a = 0; a < fr.length; a++) for (let b = a + 1; b < fr.length; b++) {
      for (const i of fr[a]) for (const j of fr[b]) {
        assert.ok(Math.hypot(X[i][0] - X[j][0], X[i][1] - X[j][1], X[i][2] - X[j][2]) > 3.0);
      }
    }
    assert.ok(/\$mol:f3/.test(r.lt));
  });

  fs.rmSync(tmp, { recursive: true, force: true });
  const failed = results.filter(r => !r[0]).length;
  console.log(`\n${results.length - failed} passed, ${failed} failed`);
  process.exit(failed ? 1 : 0);
})();
