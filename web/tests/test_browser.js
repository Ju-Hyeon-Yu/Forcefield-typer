// 빌드된 단일 HTML 을 file:// 로 열어 (네트워크 차단) 입력별 결과·복사·내려받기를 확인하는 Playwright 테스트
// 사용: node web/tests/test_browser.js [forcefield_typer.html]   (playwright 필요, 없으면 건너뜀)
const fs = require("fs");
const os = require("os");
const path = require("path");
const assert = require("assert");
const { execFileSync } = require("child_process");

let chromium;
try { ({ chromium } = require("playwright")); } catch (e) { console.log("playwright 가 없어 브라우저 테스트를 건너뜁니다"); process.exit(0); }

const HTML = path.resolve(process.argv[2] || path.join(__dirname, "../../forcefield_typer.html"));
const body = s => s.split("\n").slice(1).join("\n").replace(/^# 전하: .*$/m, "");

(async () => {
  const tmp = fs.mkdtempSync(path.join(os.tmpdir(), "dreiding-browser-"));
  execFileSync("python3", [path.join(__dirname, "make_cases.py"), tmp]);
  const cases = JSON.parse(fs.readFileSync(path.join(tmp, "cases.json")));
  const browser = await chromium.launch();
  const ctx = await browser.newContext({ viewport: { width: 1440, height: 1000 }, acceptDownloads: true });
  await ctx.grantPermissions(["clipboard-read", "clipboard-write"]);
  const blocked = [];
  await ctx.route("**/*", r => (r.request().url().startsWith("file:") ? r.continue() : (blocked.push(r.request().url()), r.abort())));
  const page = await ctx.newPage();
  const errors = [];
  page.on("pageerror", e => errors.push(e.message));
  page.on("console", m => { if (m.type() === "error") errors.push(m.text()); });
  const load = async f => {
    await page.setInputFiles("#file", f);
    await page.waitForFunction(() => !document.getElementById("fileChip").hidden);
  };
  let failed = 0;
  const test = async (name, fn) => {
    try { await fn(); console.log("ok   ", name); } catch (e) { failed++; console.log("FAIL ", name, "\n     ", e.message.split("\n")[0]); }
  };

  await test("오프라인(file://)으로 RDKit 엔진 로딩", async () => {
    await page.goto("file://" + HTML);
    await page.waitForFunction(() => document.getElementById("status").classList.contains("ready"), null, { timeout: 60000 });
  });

  await test("SMILES 만으로 .lt 생성과 타입 카드", async () => {
    await page.click("text=아세트아닐라이드");
    const lt = await page.inputValue("#lt");
    assert.ok(lt.includes('import "dreiding.lt"') && lt.includes("MOL inherits DREIDING {"));
    assert.ok(lt.includes("@atom:N_2_b1_d2_hd") && lt.includes("@atom:C_2_b1"));
    assert.ok(await page.locator(".tcard").count() >= 6);
  });

  await test(`mol2 입력 .lt 가 Python CLI 와 같음 (${cases.length}개 계)`, async () => {
    for (const c of cases) {
      await page.fill("#smiles", c.smiles);
      await load(c.mol2);
      await page.click("#run");
      assert.strictEqual(body(await page.inputValue("#lt")), body(fs.readFileSync(c.lt, "utf8")), c.smiles);
      await page.click("#clearFile");
    }
  });

  await test("GAFF2 선택: mol2 입력 .lt(보충 블록 포함)가 Python CLI 와 같음, 선택을 기억", async () => {
    await page.locator(".ffsel label", { hasText: "GAFF2" }).click();
    for (const c of cases.filter(x => x.lt_gaff2)) {
      await page.fill("#smiles", c.smiles);
      await load(c.mol2);
      await page.click("#run");
      assert.strictEqual(body(await page.inputValue("#lt")), body(fs.readFileSync(c.lt_gaff2, "utf8")), c.smiles);
      await page.click("#clearFile");
    }
    await page.fill("#smiles", "CC(N)=O");
    await page.click("#run");
    const lt = await page.inputValue("#lt");
    assert.ok(lt.includes('import "gaff2.lt"') && lt.includes("\nGAFF2 {\n") && lt.includes("MOL inherits GAFF2 {"));
    assert.ok((await page.locator(".tcard", { hasText: "nt" }).getAttribute("title")).includes("(질량 14.01)"));
    await page.reload();
    await page.waitForFunction(() => document.getElementById("status").classList.contains("ready"), null, { timeout: 60000 });
    assert.ok(await page.isChecked("input[name=ff][value=gaff2]"), "새로 열어도 GAFF2");
    await page.fill("#smiles", "F[B-](F)(F)F");
    await page.click("#run");
    assert.ok((await page.locator(".msg.err").textContent()).includes("gaff2.lt 에 없는 GAFF2 타입"));
    await page.locator(".ffsel label", { hasText: "DREIDING" }).click();
    await page.fill("#smiles", "CC(=O)Nc1ccccc1");
    await page.click("#run");
    assert.ok((await page.inputValue("#lt")).includes('import "dreiding.lt"'));
  });

  await test("OPLS-AA 선택: mol2 입력 .lt 가 Python CLI 와 같고 타입 전하를 씀, 규칙 없는 이온은 오류", async () => {
    await page.locator(".ffsel label", { hasText: "OPLS-AA" }).click();
    for (const c of cases.filter(x => x.lt_oplsaa)) {
      await page.fill("#smiles", c.smiles);
      await load(c.mol2);
      await page.click("#run");
      assert.strictEqual(body(await page.inputValue("#lt")), body(fs.readFileSync(c.lt_oplsaa, "utf8")), c.smiles);
      await page.click("#clearFile");
    }
    await page.fill("#smiles", "CC(C)(C)O");
    await page.click("#run");
    const lt = await page.inputValue("#lt");
    assert.ok(lt.includes('import "oplsaa2024.lt"') && lt.includes("MOL inherits OPLSAA {"));
    assert.ok(/@atom:160\s+-0\.650000/.test(lt), "3차 알코올 O 는 OPLS/2020 짝 160");
    assert.ok(!(await page.locator(".msg.warn").allTextContents()).some(t => t.includes("전하를 모두 0")));
    await page.fill("#smiles", "Cn1cc[n+](C)c1");
    await page.click("#run");
    assert.ok((await page.locator(".msg.err").textContent()).includes("타입이 없는 원자"));
    await page.locator(".ffsel label", { hasText: "DREIDING" }).click();
  });

  await test("pdb·xyz 입력과 전하 경고", async () => {
    for (const key of ["pdb", "pdb_noconect", "xyz"]) {
      await page.fill("#smiles", cases[0].smiles);
      await load(cases[0][key]);
      await page.click("#run");
      assert.ok((await page.locator(".msg.warn").allTextContents()).some(t => t.includes("전하 정보가 없어")), key);
      await page.click("#clearFile");
    }
  });

  await test("구조가 다르면 오류 표시", async () => {
    await page.fill("#smiles", "CCO");
    await load(cases[0].mol2);
    await page.click("#run");
    assert.ok((await page.locator(".msg.err").textContent()).includes("원자/결합 수가 다릅니다"));
    await page.click("#clearFile");
  });

  await test("타입 강조·마우스 정보·복사·내려받기", async () => {
    await page.fill("#smiles", cases[0].smiles);
    await load(cases[0].mol2);
    await page.click("#run");
    await page.locator(".tcard", { hasText: "N_R_d2" }).click();
    assert.strictEqual(await page.locator(".tcard[aria-pressed=true]").count(), 1);
    await page.click("#clearSel");
    const box = await page.locator("#view").boundingBox();
    let shown = false;
    for (let gx = 0.1; gx <= 0.9 && !shown; gx += 0.02) {
      for (let gy = 0.2; gy <= 0.8 && !shown; gy += 0.04) {
        await page.mouse.move(box.x + box.width * gx, box.y + box.height * gy);
        shown = await page.locator("#tip").isVisible();
      }
    }
    assert.ok(shown && /원자 \d+/.test(await page.textContent("#tip")), "마우스 정보");
    await page.click("#copy");
    assert.strictEqual(await page.evaluate(() => navigator.clipboard.readText()), await page.inputValue("#lt"));
    const [dl] = await Promise.all([page.waitForEvent("download"), page.click("#download")]);
    assert.strictEqual(dl.suggestedFilename(), "single_molecule.lt");
    assert.strictEqual(fs.readFileSync(await dl.path(), "utf8"), await page.inputValue("#lt"));
  });

  await test("타입 설명: 카드·원자 목록·이웃·3D 보기에서 열고, 역장을 바꿔도 같은 원자를 설명, 닫기", async () => {
    const head = () => page.textContent("#typeDetailIn h2");
    const panel = () => page.textContent("#typeDetailIn");
    await page.locator(".ffsel label", { hasText: "DREIDING" }).click();
    if (await page.locator("#fileChip").isVisible()) await page.click("#clearFile");
    await page.fill("#smiles", "CC(=O)Nc1ccccc1");
    await page.click("#run");
    assert.ok(await page.locator("#typeDetail").isHidden(), "처음에는 닫혀 있음");
    await page.click('.tcard[data-type="N_2_b1_d2_hd"]');
    assert.strictEqual(await head(), "N_2_b1_d2_hd");
    let t = await panel();
    assert.ok(t.includes("_hd: 수소결합 주개") && t.includes("판정 근거") && t.includes("아마이드·아닐린형"), t);
    assert.strictEqual(await page.locator(".main.inspect").count(), 1);
    assert.strictEqual(await page.locator(".tcard.current").count(), 1);
    // 원자 목록: 이웃 칸, 줄을 누르면 그 원자
    if (!(await page.$eval("#atomDetails", d => d.open))) await page.click("#atomSummary");
    assert.strictEqual(await page.textContent("#atomRows tr:nth-child(2) td.nbs"), "–C1 C_3, =O3 O_2_ha, –N4 N_2_b1_d2_hd");
    await page.click("#atomRows tr:nth-child(3)");
    assert.strictEqual(await head(), "O_2_ha");
    assert.strictEqual(await page.textContent(".td-atom.picked .td-aname"), "O3");
    assert.strictEqual(await page.locator("#atomRows tr.picked").count(), 1);
    // 이웃 단추로 넘어가기
    await page.locator(".td-atom.picked .nb", { hasText: "C2" }).click();
    assert.strictEqual(await head(), "C_2_b1");
    assert.ok((await panel()).includes("C2=O3 이중 ✓"));
    // 역장을 바꾸면 같은 원자(C2)의 새 타입
    await page.locator(".ffsel label", { hasText: "GAFF2" }).click();
    assert.strictEqual(await head(), "c");
    t = await panel();
    assert.ok(t.includes("ATD c * 6 3 * * [1DB,0DL] (XA1)") && t.includes("카보닐(C=O)의 sp2 C") && t.includes("Sp2 C carbonyl group"), t);
    assert.strictEqual(await page.textContent(".td-atom.picked .td-aname"), "C2");
    await page.locator(".ffsel label", { hasText: "OPLS-AA" }).click();
    assert.strictEqual(await head(), "235");
    t = await panel();
    assert.ok(t.includes("[C;X3]([O;X1])[N;X3]") && t.includes("└ [2] O, 이웃 1개"), t);
    // 3D 보기에서 원자 누르기
    const box = await page.locator("#view").boundingBox();
    let hit = null;
    for (let gy = 0.15; gy <= 0.85 && !hit; gy += 0.03) {
      for (let gx = 0.05; gx <= 0.95 && !hit; gx += 0.015) {
        const x = box.x + box.width * gx, y = box.y + box.height * gy;
        await page.mouse.move(x, y);
        if (await page.locator("#tip").isVisible()) hit = [x, y, await page.textContent("#tip strong")];
      }
    }
    assert.ok(hit, "원자를 찾지 못함");
    assert.ok((await page.textContent("#tip")).includes("눌러서 타입 설명 보기"));
    await page.mouse.click(hit[0], hit[1]);
    assert.strictEqual(await head(), hit[2]);
    assert.strictEqual(await page.locator(".td-atom.picked").count(), 1);
    // 닫기
    await page.click(".td-close");
    assert.ok(await page.locator("#typeDetail").isHidden());
    assert.strictEqual(await page.locator(".main.inspect").count(), 0);
    // 모바일 폭 확인용으로 긴 설명(패턴 나무)을 열어 둠
    await page.fill("#smiles", "c1ccncc1");
    await page.click("#run");
    await page.click('.tcard[data-type="520"]');
    assert.ok((await panel()).includes("고리 닫기: [1] 과 [6] 도 결합"));
    await page.locator(".ffsel label", { hasText: "DREIDING" }).click();
  });

  await test("모바일 폭에서 가로 스크롤 없음", async () => {
    await page.setViewportSize({ width: 390, height: 900 });
    assert.ok(await page.evaluate(() => document.documentElement.scrollWidth) <= 390);
  });

  await test("외부 네트워크 요청 없음, 페이지 오류 없음", async () => {
    assert.deepStrictEqual(blocked, []);
    assert.deepStrictEqual(errors, []);
  });

  await browser.close();
  fs.rmSync(tmp, { recursive: true, force: true });
  console.log(failed ? `\n${failed} failed` : "\nall passed");
  process.exit(failed ? 1 : 0);
})();
