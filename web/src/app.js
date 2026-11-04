// 화면 동작 — RDKit.js 준비, 입력 처리, 결과(경고·타입 카드·원자 표·뷰어·.lt) 표시
(function () {
  "use strict";
  const $ = id => document.getElementById(id);
  const T = window.DreidingTyper, P = window.DreidingPipeline;
  const tables = JSON.parse($("dreiding-tables").textContent);
  const ffs = {
    dreiding: new T.DreidingFF(tables),
    gaff2: new window.Gaff2Typer.Gaff2(JSON.parse($("gaff2-tables").textContent)),
    oplsaa: new window.FoyerTyper.Opls(JSON.parse($("oplsaa-tables").textContent)),
  };
  const HINT = {
    dreiding: "moltemplate dreiding.lt 의 라벨 규칙(결합 차수·이면각 플래그)에 맞춥니다. 전하는 mol2 의 값을 씁니다.",
    gaff2: "AmberTools antechamber(-at gaff2) 와 같은 타입을 붙입니다. gaff2.lt 에 없는 결합 파라미터는 " +
      "parmchk2 대응 타입 값으로 .lt 에 함께 씁니다. 전하는 mol2 의 값을 씁니다.",
    oplsaa: "foyer 의 OPLS-AA 규칙으로 판정해 moltemplate oplsaa2024.lt 번호를 붙입니다. 전하는 좌표 파일과 상관없이 " +
      "oplsaa2024.lt 의 타입 전하를 씁니다. 규칙이 없는 원자(이온성 액체 이온 등)는 오류로 알립니다.",
  };
  const E = window.TypeExplain;
  const viewer = new window.MolViewer($("view"), $("tip"));
  let R = null, file = null, result = null;
  let detail = { label: null, atom: -1 }, runKey = "";   // 설명 중인 타입·원자, 직전 실행의 입력 (SMILES+파일)
  viewer.onPick = i => pickAtom(i, true);
  viewer.describe = label => (result ? result.ff.describe(label) : "");

  $("srcInfo").textContent = `moltemplate ${tables.moltemplate || ""} — dreiding.lt ${Object.keys(tables.labels).length}개 라벨, ` +
    `gaff2.lt ${Object.keys(ffs.gaff2.rules.masses).length}개 타입, oplsaa2024.lt ${Object.keys(ffs.oplsaa.rules.masses).length}개 타입`;

  // ------------------------------------------------------------ 역장 선택 (브라우저에 기억)
  const ffRadios = [...document.querySelectorAll("input[name=ff]")];
  const currentFF = () => (ffRadios.find(r => r.checked) || ffRadios[0]).value;
  try {
    const saved = localStorage.getItem("fftyper.ff");
    const r = ffRadios.find(x => x.value === saved);
    if (r) r.checked = true;
  } catch (e) { /* 저장소를 못 쓰는 환경 */ }
  const showHint = () => { $("ffHint").textContent = HINT[currentFF()]; };
  showHint();
  ffRadios.forEach(r => r.addEventListener("change", () => {
    try { localStorage.setItem("fftyper.ff", currentFF()); } catch (e) { /* 무시 */ }
    showHint();
    if (R && result) run();
  }));

  // ------------------------------------------------------------ RDKit.js 준비 (압축 wasm 풀기)
  async function loadRDKit() {
    const b64 = $("rdkit-wasm").textContent.trim();
    const bin = Uint8Array.from(atob(b64), c => c.charCodeAt(0));
    if (typeof DecompressionStream === "undefined") throw new Error("이 브라우저는 압축 해제를 지원하지 않습니다. 최신 Edge 나 Chrome 을 쓰십시오");
    const stream = new Blob([bin]).stream().pipeThrough(new DecompressionStream("gzip"));
    const wasm = new Uint8Array(await new Response(stream).arrayBuffer());
    // 이 RDKit.js 빌드는 wasmBinary 옵션을 읽지 않으므로 instantiateWasm 으로 직접 넘긴다 (file:// 에서 fetch 불가)
    return window.initRDKitModule({
      instantiateWasm(imports, receive) {
        WebAssembly.instantiate(wasm, imports).then(r => receive(r.instance, r.module));
        return {};
      },
    });
  }
  loadRDKit().then(mod => {
    R = mod;
    $("status").textContent = `준비되었습니다 (RDKit ${R.version()}).`;
    $("status").classList.add("ready");
    $("run").disabled = false;
  }).catch(e => {
    $("status").textContent = "RDKit 엔진을 불러오지 못했습니다.";
    showMsgs([], `RDKit 엔진을 불러오지 못했습니다: ${e.message}`);
  });

  // ------------------------------------------------------------ 입력
  function setFile(f) {
    if (!f) { file = null; $("fileChip").hidden = true; $("drop").hidden = false; $("file").value = ""; return; }
    f.text().then(text => {
      file = { name: f.name, text };
      $("fileName").textContent = `${f.name} (${(f.size / 1024).toFixed(1)} KB)`;
      $("fileChip").hidden = false;
      $("drop").hidden = true;
    });
  }
  $("drop").addEventListener("click", () => $("file").click());
  $("drop").addEventListener("keydown", e => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); $("file").click(); } });
  $("file").addEventListener("change", e => setFile(e.target.files[0]));
  $("clearFile").addEventListener("click", () => setFile(null));
  ["dragenter", "dragover"].forEach(ev => document.addEventListener(ev, e => { e.preventDefault(); $("drop").classList.add("over"); }));
  ["dragleave", "drop"].forEach(ev => document.addEventListener(ev, e => { e.preventDefault(); if (ev === "drop" || e.target === document.documentElement) $("drop").classList.remove("over"); }));
  document.addEventListener("drop", e => { const f = e.dataTransfer && e.dataTransfer.files[0]; if (f) setFile(f); });
  $("smiles").addEventListener("keydown", e => { if (e.key === "Enter" && !$("run").disabled) run(); });
  document.querySelectorAll("[data-example]").forEach(b => b.addEventListener("click", () => {
    $("smiles").value = b.dataset.example;
    setFile(null);
    if (R) run();
  }));
  $("molName").addEventListener("input", () => {
    const v = $("molName").value.trim();
    if (/^[A-Za-z_][A-Za-z0-9_]*$/.test(v) && $("outName").dataset.touched !== "1") $("outName").value = v === "MOL" ? "single_molecule.lt" : `${v}.lt`;
  });
  $("outName").addEventListener("input", () => { $("outName").dataset.touched = "1"; });

  // ------------------------------------------------------------ 실행
  function run() {
    if (!R) return;
    try {
      result = P.run(R, ffs, { smiles: $("smiles").value, name: $("molName").value, file, ff: currentFF() });
    } catch (e) {
      if (!(e instanceof T.TypingError) && !(e instanceof window.DreidingIO.InputError)) console.error(e);
      result = null;
      clearDetail();
      showMsgs([], e.message);
      return;
    }
    showMsgs(result.warnings, null);
    render(result);
  }
  $("run").addEventListener("click", run);

  function showMsgs(warnings, error) {
    const box = $("msgs");
    box.textContent = "";
    if (error) {
      const d = document.createElement("div");
      d.className = "msg err";
      d.textContent = error;
      box.append(d);
    }
    for (const w of warnings) {
      const d = document.createElement("div");
      d.className = "msg warn";
      d.textContent = w;
      box.append(d);
    }
  }

  function ranges(nums) {
    const out = [];
    for (let k = 0; k < nums.length; k++) {
      let j = k;
      while (j + 1 < nums.length && nums[j + 1] === nums[j] + 1) j++;
      out.push(j - k >= 2 ? `${nums[k]}–${nums[j]}` : j > k ? `${nums[k]}, ${nums[j]}` : `${nums[k]}`);
      k = j;
    }
    return out.join(", ");
  }

  function render(r) {
    // 타입 순서: 원소 기호 순, 같은 원소 안에서는 처음 나온 순
    const firstSeen = new Map();
    r.rows.forEach((row, i) => { if (!firstSeen.has(row.label)) firstSeen.set(row.label, i); });
    const types = [...firstSeen.keys()].sort((a, b) => {
      const ea = r.rows[firstSeen.get(a)].sym, eb = r.rows[firstSeen.get(b)].sym;
      const rank = e => (e === "C" ? 0 : e === "H" ? 1 : 2);
      return rank(ea) - rank(eb) || ea.localeCompare(eb) || firstSeen.get(a) - firstSeen.get(b);
    });
    const flagged = new Set(r.flagged);
    viewer.setMolecule({ atoms: r.rows, bonds: r.bonds, flagged }, types);
    $("empty").hidden = true;
    $("toolbar").hidden = false;
    $("legend").hidden = flagged.size === 0;

    const grid = $("typeGrid");
    grid.textContent = "";
    for (const t of types) {
      const idx = r.rows.map((row, i) => (row.label === t ? i + 1 : 0)).filter(Boolean);
      const b = document.createElement("button");
      b.type = "button";
      b.className = "tcard";
      b.setAttribute("aria-pressed", "false");
      b.dataset.type = t;
      const sw = document.createElement("span");
      sw.className = "sw";
      sw.style.background = viewer.typeColor.get(t);
      const nm = document.createElement("span");
      nm.className = "name";
      nm.textContent = t;
      const ct = document.createElement("span");
      ct.className = "count";
      ct.textContent = `${idx.length}개`;
      const at = document.createElement("span");
      at.className = "atoms";
      at.textContent = `원자 ${ranges(idx)}`;
      const ds = document.createElement("span");
      ds.className = "desc";
      ds.textContent = r.ff.describe(t);
      b.append(sw, nm, ct, at, ds);
      b.title = `${r.ff.describe(t)} — 눌러서 강조하고 타입 설명 보기`;
      b.addEventListener("click", () => { toggleType(t, b); showDetail(t, detail.label === t ? detail.atom : -1); });
      grid.append(b);
    }
    $("typeSummary").textContent = `원자 ${r.rows.length}개, 타입 ${types.length}종, 조각 ${r.nFrags}개입니다. ` +
      "타입이나 원자를 누르면 그 타입의 뜻과 판정 근거를 보여 줍니다.";
    $("clearSel").hidden = true;

    const tb = $("atomRows");
    tb.textContent = "";
    r.rows.forEach((row, i) => {
      const tr = document.createElement("tr");
      if (flagged.has(i)) tr.className = "flag";
      const cells = [[i + 1, "num"], [row.name, ""], [row.label, "lab"], [E.neighborText(r, i), "nbs"], [row.q.toFixed(6), "num"],
        [row.frag, "num"], [row.xyz[0].toFixed(3), "num"], [row.xyz[1].toFixed(3), "num"], [row.xyz[2].toFixed(3), "num"]];
      for (const [v, cls] of cells) { const td = document.createElement("td"); td.textContent = v; if (cls) td.className = cls; tr.append(td); }
      tr.tabIndex = 0;
      tr.title = "눌러서 이 원자의 타입 설명 보기";
      tr.addEventListener("click", () => pickAtom(i, false));
      tr.addEventListener("keydown", e => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); pickAtom(i, false); } });
      tb.append(tr);
    });
    $("atomSummary").textContent = `원자별 목록 (${r.rows.length}개)`;
    $("atomDetails").hidden = false;

    $("lt").value = r.lt;
    $("ltTitle").textContent = `${outName()} (${r.name}, ${r.ff.LABEL})`;
    $("copy").disabled = false;
    $("download").disabled = false;
    $("copied").textContent = "";

    // 같은 입력을 다른 역장으로 다시 돌린 경우: 보던 원자(또는 같은 이름의 타입) 설명을 이어서 보여 줌
    const key = `${$("smiles").value.trim()}\u0000${file ? file.name + file.text.length : ""}`;
    const prev = detail;
    clearDetail();
    if (key === runKey && prev.atom >= 0 && prev.atom < r.rows.length) showDetail(r.labels[prev.atom], prev.atom);
    else if (key === runKey && prev.label !== null && r.labels.includes(prev.label)) showDetail(prev.label);
    runKey = key;
  }

  // ------------------------------------------------------------ 타입 설명
  const wide = () => window.matchMedia("(min-width: 1100px)").matches;
  const mk = (tag, cls, text) => {
    const e = document.createElement(tag);
    if (cls) e.className = cls;
    if (text !== undefined) e.textContent = text;
    return e;
  };

  function lineList(lines) {
    const ul = mk("ul", "td-lines");
    for (const l of lines) {
      if (typeof l === "string") {
        const li = mk("li", l.includes("✗") ? "bad" : "", l);
        ul.append(li);
      } else {
        const li = mk("li", "tree");
        li.append(mk("pre", "td-tree", l.tree.join("\n")));
        ul.append(li);
      }
    }
    return ul;
  }

  function atomBlock(a, picked) {
    const box = mk("div", "td-atom" + (picked ? " picked" : ""));
    const head = mk("div", "td-ahead");
    const nm = mk("button", "td-aname", a.name);
    nm.type = "button";
    nm.title = "이 원자를 3D 보기와 원자 목록에서 표시";
    nm.addEventListener("click", () => pickAtom(a.i, false, true));
    head.append(nm);
    if (a.flagged) head.append(mk("span", "td-flag", "경고 대상"));
    box.append(head);
    const nb = mk("div", "td-nbrs");
    if (!a.env.nbrs.length) nb.textContent = "이웃 없음 (홀로 있는 원자·이온)";
    else {
      nb.append(`이웃 ${a.env.nbrs.length}개: `);
      for (const x of a.env.nbrs) {
        const b = mk("button", "nb");
        b.type = "button";
        b.title = `${x.name} (${x.label}) 의 타입 설명으로`;
        b.append(`${E.BOND_SYM[x.order] || "–"}${x.name} `, mk("b", "", x.label));
        b.addEventListener("click", () => pickAtom(x.j, false, true));
        nb.append(b);
      }
    }
    box.append(nb);
    const env = E.envText(a.env);
    if (env.length > 1) box.append(mk("div", "td-env", env[1]));
    if (a.lines.length) box.append(lineList(a.lines));
    return box;
  }

  function showDetail(label, atom = -1, focus = false) {
    if (!result) return;
    detail = { label, atom };
    const x = E.explainType(result, label);
    const box = $("typeDetailIn");
    box.textContent = "";

    const head = mk("div", "td-head");
    const sw = mk("span", "sw");
    sw.style.background = viewer.typeColor.get(label) || "#999";
    const h = mk("h2", "", label);
    h.tabIndex = -1;
    const close = mk("button", "td-close", "×");
    close.type = "button";
    close.title = "타입 설명 닫기";
    close.setAttribute("aria-label", "타입 설명 닫기");
    close.addEventListener("click", clearDetail);
    head.append(sw, h, close, mk("p", "sub", `${x.ffLabel} 타입 · 원자 ${x.atomIdx.length}개 (${ranges(x.atomIdx.map(i => i + 1))})`));
    box.append(head);

    const sec1 = mk("section", "td-sec");
    sec1.append(mk("h3", "", "뜻"));
    if (x.meaning.original) {
      const q = mk("blockquote", "td-orig");
      q.append(mk("small", "", `원문 · ${x.meaning.source}`), x.meaning.original);
      sec1.append(q);
    }
    sec1.append(lineList(x.meaning.lines));
    box.append(sec1);

    const sec2 = mk("section", "td-sec");
    sec2.append(mk("h3", "", "판정 근거"));
    for (const g of x.groups) {
      const gb = mk("div", "td-group");
      const t = mk("div", "td-gtitle", g.title + " ");
      if (x.groups.length > 1 || x.atomIdx.length > 1) t.append(mk("span", "", `· 원자 ${ranges(g.atoms.map(i => i + 1))}`));
      gb.append(t);
      if (g.code) gb.append(mk("code", "td-code", g.code));
      gb.append(lineList(g.lines));
      sec2.append(gb);
    }
    box.append(sec2);

    const sec3 = mk("section", "td-sec");
    const h3 = mk("h3", "", "원자별 실제 환경 ");
    h3.append(mk("span", "", "— 이웃을 누르면 그 원자의 설명으로 넘어갑니다"));
    sec3.append(h3);
    const order = x.atoms.slice().sort((p, q) => (q.i === atom) - (p.i === atom));
    const SHOW = 8;
    order.slice(0, SHOW).forEach(a => sec3.append(atomBlock(a, a.i === atom)));
    if (order.length > SHOW) {
      const more = mk("button", "link td-more", `나머지 원자 ${order.length - SHOW}개 더 보기`);
      more.type = "button";
      more.addEventListener("click", () => { order.slice(SHOW).forEach(a => sec3.insertBefore(atomBlock(a, false), more)); more.remove(); });
      sec3.append(more);
    }
    box.append(sec3);

    $("typeDetail").hidden = false;
    document.querySelector(".main").classList.add("inspect");
    viewer.setPicked(atom);
    document.querySelectorAll(".tcard").forEach(b => b.classList.toggle("current", b.dataset.type === label));
    document.querySelectorAll("#atomRows tr").forEach((tr, k) => tr.classList.toggle("picked", k === atom));
    box.scrollTop = 0;
    if (focus) h.focus({ preventScroll: true });
  }

  /** fromView: 3D 보기에서 누름 (좁은 화면이면 설명으로 내려감), focus: 설명 안에서 누름 (제목으로 초점 이동). */
  function pickAtom(i, fromView, focus = false) {
    if (!result) return;
    showDetail(result.labels[i], i, focus);
    if (fromView && !wide()) $("typeDetail").scrollIntoView({ block: "start", behavior: "smooth" });
    // 원자 목록이 열려 있으면 그 줄이 보이도록 표 안에서만 굴림 (페이지는 그대로)
    const tr = $("atomRows").children[i], wrap = tr && tr.closest(".tablewrap");
    if (wrap && $("atomDetails").open) {
      const head = wrap.querySelector("thead").offsetHeight, top = tr.offsetTop, bottom = top + tr.offsetHeight;
      if (top - head < wrap.scrollTop) wrap.scrollTop = top - head;
      else if (bottom > wrap.scrollTop + wrap.clientHeight) wrap.scrollTop = bottom - wrap.clientHeight;
    }
  }

  function clearDetail() {
    detail = { label: null, atom: -1 };
    $("typeDetail").hidden = true;
    $("typeDetailIn").textContent = "";
    document.querySelector(".main").classList.remove("inspect");
    viewer.setPicked(-1);
    document.querySelectorAll(".tcard.current").forEach(b => b.classList.remove("current"));
    document.querySelectorAll("#atomRows tr.picked").forEach(tr => tr.classList.remove("picked"));
  }

  function toggleType(t, btn) {
    if (viewer.selected.has(t)) viewer.selected.delete(t); else viewer.selected.add(t);
    btn.setAttribute("aria-pressed", String(viewer.selected.has(t)));
    $("clearSel").hidden = viewer.selected.size === 0;
    viewer.draw();
  }
  $("clearSel").addEventListener("click", () => {
    viewer.selected.clear();
    document.querySelectorAll(".tcard").forEach(b => b.setAttribute("aria-pressed", "false"));
    $("clearSel").hidden = true;
    viewer.draw();
  });

  // ------------------------------------------------------------ 뷰어 도구
  document.querySelectorAll("[data-label]").forEach(b => b.addEventListener("click", () => {
    viewer.labelMode = b.dataset.label;
    document.querySelectorAll("[data-label]").forEach(x => x.setAttribute("aria-pressed", String(x === b)));
    viewer.draw();
  }));
  document.querySelectorAll("[data-color]").forEach(b => b.addEventListener("click", () => {
    viewer.colorMode = b.dataset.color;
    document.querySelectorAll("[data-color]").forEach(x => x.setAttribute("aria-pressed", String(x === b)));
    viewer.draw();
  }));
  $("resetView").addEventListener("click", () => viewer.resetView());
  $("png").addEventListener("click", () => {
    const src = $("view"), c = document.createElement("canvas");
    c.width = src.width; c.height = src.height;
    const ctx = c.getContext("2d");
    ctx.fillStyle = getComputedStyle(document.querySelector(".stage")).backgroundColor;
    ctx.fillRect(0, 0, c.width, c.height);
    ctx.drawImage(src, 0, 0);
    c.toBlob(b => save(b, outName().replace(/\.lt$/, "") + ".png"));
  });

  // ------------------------------------------------------------ 출력
  function outName() {
    let n = $("outName").value.trim() || "single_molecule.lt";
    if (!/\.lt$/i.test(n)) n += ".lt";
    return n.replace(/[\\/:*?"<>|]/g, "_");
  }
  function save(blob, name) {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = name;
    document.body.append(a);
    a.click();
    setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
  }
  $("download").addEventListener("click", () => {
    if (result) save(new Blob([result.lt], { type: "text/plain;charset=utf-8" }), outName());
  });
  $("copy").addEventListener("click", async () => {
    if (!result) return;
    let ok = false;
    try { await navigator.clipboard.writeText(result.lt); ok = true; } catch (e) {
      const ta = $("lt");
      ta.focus(); ta.select();
      try { ok = document.execCommand("copy"); } catch (e2) { ok = false; }
    }
    $("copied").textContent = ok ? "복사했습니다." : "복사하지 못했습니다. 아래 칸을 직접 선택해 복사하십시오.";
  });
})();
