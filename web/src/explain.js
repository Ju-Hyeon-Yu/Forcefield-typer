// 할당된 원자 타입의 뜻·판정 근거·실제 원자 환경을 화면에 보여 줄 문장으로 만드는 모듈 (DREIDING·GAFF2·OPLS-AA)
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory(require("./typer.js"));
  else root.TypeExplain = factory(root.DreidingTyper);
})(typeof self !== "undefined" ? self : this, function (T) {
  "use strict";

  /*
   * explainType(result, label) → {
   *   label, ffLabel, atomIdx: [i],
   *   meaning: { original, source, lines },   // original: 역장 파일의 원문 설명(영어), lines: 한국어 풀이
   *   groups: [{ title, code, lines, atoms }], // 판정 근거 — 같은 근거로 타입을 받은 원자끼리 묶음
   *   atoms: [{ i, name, flagged, env, lines }] // 원자별 실제 환경 (env 는 이웃·H·고리, lines 는 역장별 판정 값)
   * }
   * result 는 pipeline.run 의 반환값 (graph·labels·ff·info·flagged).
   */

  const BOND_SYM = { 1: "–", 2: "=", 3: "≡", 1.5: ":" };
  const ORDER_KO = { 1: "단일", 1.5: "방향족(1.5)", 2: "이중", 3: "삼중" };
  const atomName = (g, i) => `${g.atoms[i].sym}${i + 1}`;

  // ---------------------------------------------------------------- 공통: 실제 원자 환경
  /** 원자마다 그 원자를 지나는 가장 작은 고리의 크기 (고리가 아니면 0). */
  function smallestRings(g) {
    const n = g.atoms.length, out = new Array(n).fill(0);
    for (let i = 0; i < n; i++) {
      const nb = g.nbrs[i].map(([j]) => j);
      let best = 0;
      for (let p = 0; p < nb.length; p++) {
        // i 를 지나지 않고 nb[p] 에서 다른 이웃까지 가는 가장 짧은 길 + 2 = 고리 크기
        const dist = new Map([[nb[p], 0]]), queue = [nb[p]];
        for (let h = 0; h < queue.length; h++) {
          const u = queue[h];
          for (const [v] of g.nbrs[u]) if (v !== i && !dist.has(v)) { dist.set(v, dist.get(u) + 1); queue.push(v); }
        }
        for (let q = p + 1; q < nb.length; q++) {
          if (!dist.has(nb[q])) continue;
          const size = dist.get(nb[q]) + 2;
          if (!best || size < best) best = size;
        }
      }
      out[i] = best;
    }
    return out;
  }

  function environment(r, i) {
    const g = r.graph;
    if (!r._ringSize) r._ringSize = smallestRings(g);
    return {
      sym: g.atoms[i].sym,
      nbrs: g.nbrs[i].map(([j, o]) => ({ j, name: atomName(g, j), label: r.labels[j], order: o })),
      nh: g.nbrs[i].filter(([j]) => g.atoms[j].sym === "H").length,
      chg: g.atoms[i].chg,
      ring: r._ringSize[i],
      aromatic: g.atoms[i].aromatic,
    };
  }

  /** 환경을 글로 (테스트·복사용). */
  function envText(env) {
    const nb = env.nbrs.map(x => `${BOND_SYM[x.order] || "–"}${x.name} (${x.label})`);
    const head = nb.length ? `이웃 ${nb.length}개: ${nb.join(", ")}` : "이웃 없음 (홀로 있는 원자·이온)";
    if (env.sym === "H") return [head];
    const parts = [`H ${env.nh}개`];
    if (env.chg) parts.push(`형식 전하 ${env.chg > 0 ? "+" : ""}${env.chg}`);
    parts.push(env.ring ? `고리 원자 (가장 작은 고리 ${env.ring}원${env.aromatic ? ", 방향족" : ""})` : "고리 아님");
    return [head, parts.join(" · ")];
  }

  /** 원자 목록 표의 "이웃 (타입)" 칸. */
  function neighborText(r, i) {
    const g = r.graph;
    return g.nbrs[i].map(([j, o]) => `${BOND_SYM[o] || "–"}${atomName(g, j)} ${r.labels[j]}`).join(", ");
  }

  // ---------------------------------------------------------------- DREIDING
  const HYB = { 3: "sp3 (사면체)", 2: "sp2 (평면 삼각형)", 1: "sp (직선)", R: "공명 — 방향족 고리 원자 (평면)" };
  const FLAG_RE = /^([A-Z][a-z]?_[123R]\d?)(_b[12])?(.*)$/;

  function dreidingMeaning(label, ff) {
    if (label === "H") return ["일반 H — 수소결합을 하지 않는 H (C 등에 붙음)"];
    if (label === "H_HB") return ["수소결합 H — N·O·F 에 붙어 수소결합 주개가 되는 H"];
    if (label === "H_B") return ["붕소 화합물의 다리 H"];
    const p = label.split("_");
    const lines = [`원소 ${p[0]}`];
    for (const part of p.slice(1)) {
      let m;
      if ((m = /^([123R])(\d?)$/.exec(part))) {
        lines.push(`_${part}: ${HYB[m[1]]}` + (m[2] ? `, H ${m[2]}개를 품은 united-atom 라벨` : ""));
      } else if ((m = /^b([12])$/.exec(part))) {
        lines.push(`_b${m[1]}: 결합 차수 플래그 — sp2 원자끼리의 결합 차수를 두 라벨의 플래그 조합으로 정함 (반드시 짝으로 씀)`);
      } else if ((m = /^d(\d)$/.exec(part))) {
        lines.push(`_d${m[1]}: 이웃 ${Number(m[1]) + 1}개용 라벨`);
      } else if (part === "hd") {
        lines.push("_hd: 수소결합 주개 (H 가 붙은 N·O·F)");
      } else if (part === "ha") {
        lines.push("_ha: 수소결합 받개");
      }
    }
    if (p.length === 1 && ["F", "Cl", "Br", "I"].includes(p[0])) lines.push("할로젠 (이웃 1개용)");
    else if (p.length === 1) lines.push("홀로 있는 이온");
    const fm = FLAG_RE.exec(label);
    if (fm && T.sp2Class(label)) {
      // 플래그의 효과: 같은 혼성의 C 라벨(플래그 없음·_b1·_b2)과 결합할 때의 차수
      const hyb = fm[1].split("_")[1][0];
      const pairs = ["", "_b1", "_b2"].map(f => `C_${hyb}${f}`).filter(l => l in ff.labels).map(v => {
        const o = ff.bondOrder(label, v);
        return `${v} → ${o === null ? "규칙 없음" : ORDER_KO[o] || o}`;
      });
      if (pairs.length) lines.push(`C_${hyb} 계열과 결합할 때 moltemplate 가 쓰는 차수: ${pairs.join(", ")}`);
    }
    const t0 = ff.theta0 && ff.theta0[label];
    if (t0 !== undefined) lines.push(`평형 결합각 θ0 = ${t0}° (dreiding.lt)`);
    return lines;
  }

  /** typer.js baseLabel 의 갈림길을 같은 순서로 따라가며 이유를 문장으로. */
  function dreidingBaseReason(g, i) {
    const a = g.atoms[i], sym = a.sym, deg = g.deg[i];
    const orders = g.nbrs[i].map(([, o]) => o);
    const has2 = orders.includes(2), has3 = orders.includes(3);
    if (sym === "H") {
      const h = g.nbrs[i][0] && g.atoms[g.nbrs[i][0][0]].sym;
      return ["N", "O", "F"].includes(h) && deg === 1 ? `${h} 에 붙은 H → H_HB` : "N·O·F 가 아닌 원자에 붙은 H → H";
    }
    if (sym === "C") {
      if (a.aromatic) return "방향족 고리 원자 → C_R";
      return { 4: "이웃 4개 → sp3, C_3", 3: "이웃 3개 → sp2, C_2", 2: "이웃 2개 → sp, C_1" }[deg] || `이웃 ${deg}개`;
    }
    if (sym === "N") {
      if (a.aromatic) return `방향족 고리 N → N_R, 이웃 ${deg}개 → _d${deg - 1}`;
      if (deg === 4) return "이웃 4개 (암모늄형) → N_3";
      if (has3) return "삼중결합 N → N_1";
      if (has2) return `이중결합 N → N_2, 이웃 ${deg}개 → _d${deg - 1}`;
      if (g.nbrs[i].some(([n]) => T.isSp2Atom(g, n))) {
        return `C=O·C=C·방향족 고리 같은 sp2 원자에 붙은 N (아마이드·아닐린형, 공액) → N_2, 이웃 ${deg}개 → _d${deg - 1}`;
      }
      return "sp2 이웃이 없는 N (아민형) → N_3";
    }
    if (sym === "O") {
      if (a.aromatic) return "방향족 고리 O → O_R";
      if (deg === 1) {
        const n = g.nbrs[i][0][0], nb = g.atoms[n].sym;
        if (nb === "S" || nb === "P") return `${nb} 에 붙은 말단 O → O_2`;
        if (has2) return "이중결합 말단 O (C=O 등) → O_2";
        if (has3) return "삼중결합 O → O_1";
        if ((nb === "C" || nb === "N") && !g.atoms[n].aromatic && T.isSp2Atom(g, n)) {
          return "방향족이 아닌 sp2 원자에 붙은 말단 O (카복실레이트·니트로의 –O⁻) → O_2";
        }
        return g.atoms[n].aromatic ? "방향족 고리 원자에 붙은 말단 O (페놀레이트형) → O_3" : "sp3 원자에 붙은 말단 O (알콕사이드형) → O_3";
      }
      if (deg === 3) return "3배위 O — 맞는 라벨이 없어 O_3 (경고 대상)";
      return "이웃 2개 (알코올·에터·에스터 O) → O_3";
    }
    if (sym === "B") return deg === 4 ? "이웃 4개 → B_3" : `이웃 ${deg}개 → B_2, _d${deg - 1}`;
    if (["F", "Cl", "Br", "I"].includes(sym)) return `할로젠 → ${sym}`;
    if (sym === "S") return `dreiding.lt 의 S 라벨은 S_3 뿐 (이 S 는 이웃 ${deg}개)`;
    if (sym === "P") return `이웃 ${deg}개 → P_3_d${deg - 1}`;
    if (sym === "Si") return "Si → Si_3";
    return deg === 0 ? `홀로 있는 이온 → ${sym}` : `이웃 ${deg}개`;
  }

  /** typer.js hFlag 와 같은 판단의 이유. */
  function dreidingHReason(g, i) {
    const a = g.atoms[i], sym = a.sym;
    if (sym !== "N" && sym !== "O" && sym !== "F") return null;
    if (g.nbrs[i].some(([n]) => g.atoms[n].sym === "H")) return "H 가 붙어 있음 → _hd (수소결합 주개)";
    if (sym === "O" || sym === "F") return `H 가 없는 ${sym} → _ha (수소결합 받개)`;
    const deg = g.deg[i];
    if (a.chg > 0) return "양전하 N → 수소결합 플래그 없음";
    if (a.aromatic) return deg === 2 ? "이웃 2개인 방향족 N (피리딘형) → _ha" : "이웃 3개인 방향족 N (고립 전자쌍이 고리에 쓰임) → 플래그 없음";
    const orders = g.nbrs[i].map(([, o]) => o);
    if (orders.includes(3) || orders.includes(2) || deg === 2) return "다중결합 N 이나 이웃 2개인 N → _ha";
    if (deg === 3 && !g.nbrs[i].some(([n]) => T.isSp2Atom(g, n))) return "아민 N → _ha";
    return "sp2 원자에 붙은 N (고립 전자쌍이 공액에 쓰임) → 수소결합 플래그 없음";
  }

  function dreidingExplain(r, label, idx) {
    const g0 = r.graph, ff = r.ff;
    const [g, dearo] = T.dearomatizeExocyclic(g0);
    const dset = new Set(dearo);
    const groups = new Map();
    for (const i of idx) {
      const lines = [];
      if (dset.has(i)) lines.push("고리 밖 이중결합(C=O 등)이 붙은 방향족 고리라 비방향족 sp2 고리(_2 + b 플래그)로 다룸");
      lines.push(dreidingBaseReason(g, i));
      const hr = dreidingHReason(g, i);
      if (hr) lines.push(hr);
      const want = T.hFlag(g, i);
      if (want && !label.endsWith(want)) lines.push(`'${label}${want}' 라벨이 dreiding.lt 에 없어 ${want} 를 뺌 (경고 대상)`);
      if (/_b[12]/.test(label)) {
        lines.push("_b 플래그는 결합 차수가 SMILES 와 같아지고, 모든 결합에 이면각 규칙이 있고, 같은 플래그 이웃이 있도록 " +
          "고른 조합 (플래그가 적은 쪽을 먼저 고름)");
      }
      const key = lines.join("\n");
      if (!groups.has(key)) groups.set(key, { title: "라벨을 고른 이유", code: null, lines, atoms: [] });
      groups.get(key).atoms.push(i);
    }
    const atoms = idx.map(i => {
      const lines = [];
      if (T.sp2Class(label)) {
        // 기대 차수는 고리 밖 이중결합 처리를 한 그래프 기준 (타이퍼와 같음)
        const chk = g.nbrs[i].filter(([j]) => g.atoms[j].sym !== "H").map(([j, want]) => {
          const got = ff.bondOrder(r.labels[i], r.labels[j]);
          const ok = got !== null && Math.abs(got - want) < 1e-6;
          return `${atomName(g, i)}${BOND_SYM[want] || "–"}${atomName(g, j)} ${ORDER_KO[want] || want}` +
            (ok ? " ✓" : ` → ${got === null ? "규칙 없음" : ORDER_KO[got] || got} ✗`);
        });
        if (chk.length) lines.push("결합 차수 (SMILES 기준 → moltemplate 라벨 쌍): " + chk.join(", "));
      }
      return { i, name: atomName(g0, i), flagged: r.flagged.includes(i), env: environment(r, i), lines };
    });
    return { meaning: { original: null, source: null, lines: dreidingMeaning(label, ff) }, groups: [...groups.values()], atoms };
  }

  // ---------------------------------------------------------------- GAFF2
  // gaff2.lt 원문 설명의 한국어 풀이
  const GAFF2_KO = {
    c: "카보닐(C=O)의 sp2 C", cs: "C=S 의 sp2 C", c1: "sp C", c2: "sp2 C (일반)", c3: "sp3 C",
    ca: "순수 방향족 고리(벤젠·피리딘형)의 sp2 C", cp: "두 방향족 고리를 잇는 바이페닐형 연결 sp2 C",
    cq: "cp 와 같은 원자 — 짝 표기용 이름", cc: "순수 방향족이 아닌 고리 공액계의 sp2 C",
    cd: "cc 와 같은 원자 — 공액 짝 표기용 이름", ce: "고리 밖 공액계 안쪽의 sp2 C", cf: "ce 와 같은 원자 — 공액 짝 표기용 이름",
    cg: "공액계 안쪽의 sp C", ch: "cg 와 같은 원자 — 공액 짝 표기용 이름", cx: "3원 고리의 sp3 C", cy: "4원 고리의 sp3 C",
    cu: "3원 고리의 sp2 C", cv: "4원 고리의 sp2 C", cz: "구아니딘기의 sp2 C",
    h1: "전자 끄는 원자 1개가 붙은 지방족 C 의 H", h2: "전자 끄는 원자 2개가 붙은 지방족 C 의 H",
    h3: "전자 끄는 원자 3개가 붙은 지방족 C 의 H", h4: "전자 끄는 원자 1개가 붙은 sp3 이 아닌 C 의 H",
    h5: "전자 끄는 원자 2개가 붙은 sp3 이 아닌 C 의 H", ha: "방향족 C 의 H", hc: "전자 끄는 원자가 없는 지방족 C 의 H",
    hn: "N 에 붙은 H", ho: "하이드록실 H", hp: "P 에 붙은 H", hs: "S 에 붙은 H", hw: "물의 H", hx: "양전하 기 옆 C 의 H",
    f: "플루오린", cl: "염소", br: "브로민", i: "아이오딘",
    n: "아마이드기의 sp2 N", n1: "sp N", n2: "이웃 2개인 지방족 sp2 N", n3: "이웃 3개인 sp3 N", n4: "이웃 4개인 sp3 N",
    na: "이웃 3개인 sp2 N", nb: "순수 방향족 고리(피리딘형)의 sp2 N", nc: "순수 방향족이 아닌 고리의 sp2 N",
    nd: "nc 와 같은 원자 — 공액 짝 표기용 이름", ne: "공액계 안쪽의 sp2 N", nf: "ne 와 같은 원자 — 공액 짝 표기용 이름",
    nh: "방향족 고리에 붙은 아민 N", no: "니트로 N", ns: "H 1개가 붙은 아마이드 N", nt: "H 2개가 붙은 아마이드 N",
    nx: "n4 와 같으나 H 1개", ny: "n4 와 같으나 H 2개", nz: "n4 와 같으나 H 3개", "n+": "NH4+ 의 N",
    nu: "nh 와 같으나 H 1개", nv: "nh 와 같으나 H 2개", n7: "n3 와 같으나 H 1개", n8: "n3 와 같으나 H 2개", n9: "NH3 의 N",
    ni: "3원 고리의 n", nj: "4원 고리의 n", nk: "3원 고리의 n4·nx·ny", nl: "4원 고리의 n4·nx·ny", nm: "3원 고리의 nh",
    nn: "4원 고리의 nh", np: "3원 고리의 n3", nq: "4원 고리의 n3", n5: "3원 고리의 n7", n6: "4원 고리의 n7",
    o: "이웃 1개인 O (카보닐·카복실레이트 등)", oh: "하이드록실 O", os: "에터·에스터 O", op: "3원 고리의 os", oq: "4원 고리의 os",
    ow: "물의 O", p2: "이웃 2개인 P", p3: "이웃 3개인 P (PH3 형)", p4: "이웃 3개인 P (O=P(CH3)2 형)",
    p5: "이웃 4개인 P (O=P(OH)3 형)", pb: "순수 방향족 고리의 sp2 P", pc: "순수 방향족이 아닌 고리의 sp2 P",
    pd: "pc 와 같은 원자 — 공액 짝 표기용 이름", pe: "공액계 안쪽의 sp2 P", pf: "pe 와 같은 원자 — 공액 짝 표기용 이름",
    px: "공액계의 특수 p4", py: "공액계의 특수 p5", s: "이웃 1개인 S", s2: "이웃 2개이고 이중결합이 있는 S",
    s4: "이웃 3개인 S", s6: "이웃 4개인 S", sh: "H 가 붙은 sp3 S", ss: "싸이오에터·싸이오에스터의 sp3 S",
    sp: "3원 고리의 ss", sq: "4원 고리의 ss", sx: "공액계의 특수 s4", sy: "공액계의 특수 s6", Cu: "구리",
  };

  // ATOMTYPE_GFF2.DEF 의 미리 정한 낱말 (파일 끝 설명을 옮김)
  const AP_KO = {
    RG: "고리", NR: "고리 밖", AR1: "AR1 고리(순수 방향족)", AR2: "AR2 고리(평면, 단일결합 둘이 잇달고 이중결합 둘 이상)",
    AR3: "AR3 고리(평면, 고리 밖 원자와 이중결합)", AR4: "AR4 고리(그 밖)", AR5: "AR5 고리(sp3 C 로만 된 지방족)",
    SB: "단일결합(방향족 제외)", sb: "단일결합(방향족 단일·비편재 포함)", DB: "이중결합(방향족 제외)",
    db: "이중결합(방향족 이중 포함)",
    TB: "삼중결합", tb: "삼중결합", AB: "방향족 결합", DL: "비편재 결합",
  };
  const WILD_KO = { XX: "C·N·O·S·P", XA: "O·S", XB: "N·P", XC: "할로젠", XD: "S·P" };
  const GBOND_KO = { 1: "단일", 2: "이중", 3: "삼중", 6: "N–O 말단(결합 셈에서 빠짐)", 7: "방향족 단일", 8: "방향족 이중",
    9: "비편재", 10: "방향족" };

  function apUnitKo([count, name]) {
    let m = /^RG(\d+)$/.exec(name);
    if (m) return count === -1 ? `${m[1]}원 고리에 속함` : count === 0 ? `${m[1]}원 고리에 속하지 않음` : `${m[1]}원 고리 ${count}개에 속함`;
    if (name === "RG") return count === -1 ? "고리 원자" : count === 0 ? "고리 원자 아님" : `고리 ${count}개에 속함`;
    if (name === "NR") return count === 0 ? "고리 원자" : "고리 밖 원자";
    m = /^AR([1-5])$/.exec(name);
    if (m) return count === -1 ? `${AP_KO[name]}에 속함` : count === 0 ? `${AP_KO[name]}에 속하지 않음` : `${AP_KO[name]} ${count}개에 속함`;
    m = /^([A-Za-z]{2})('*)$/.exec(name);
    if (m && AP_KO[m[1]]) {
      const base = AP_KO[m[1]];
      let s = count === -1 ? `${base} 있음` : count === 0 ? `${base} 없음` : `${base} ${count}개`;
      if (m[2] === "'") s += " — 바로 앞 원자와의 결합이 이 종류";
      else if (m[2] === "''") s += " — 바로 앞 원자와의 결합은 이 종류가 아님";
      return s;
    }
    return `${count === -1 ? "" : count}${name}`;
  }
  const apKo = ap => ap.map(alts => alts.map(apUnitKo).join(" 또는 ")).join(", ");
  const parseApText = text => text.replace(/^[\[\]]+|[\[\]]+$/g, "").split(",").map(c => c.split(".").map(u => {
    const m = /^(\d*)([\s\S]*)$/.exec(u);
    return [m[1] ? parseInt(m[1], 10) : -1, m[2]];
  }));

  function cesKo(chains) {
    return chains.map((ch, k) => `이웃 경로 ${k + 1}: 이 원자 → ` + ch.map(nd => {
      const el = WILD_KO[nd.name] ? `${nd.name}(${WILD_KO[nd.name]} 중 하나)` : nd.name;
      const con = nd.connum ? `, 이웃 ${nd.connum}개` : "";
      const ap = nd.ap ? ` [${apKo(parseApText(nd.ap))}]` : "";
      return `${el}${con}${ap}`;
    }).join(" → "));
  }

  function gaffRuleLines(rule) {
    const f = [];
    if (rule.z !== undefined && rule.z !== null) f.push(`원자번호 ${rule.z} (${T.SYMBOLS[rule.z]})`);
    if (rule.connum !== undefined && rule.connum !== null) f.push(`이웃 ${rule.connum}개`);
    if (rule.nh !== undefined && rule.nh !== null) f.push(`붙은 H ${rule.nh}개`);
    if (rule.new !== undefined && rule.new !== null) f.push(`붙은 원자에 달린 전자 끄는 원자(N·O·S·할로젠) ${rule.new}개`);
    const out = [f.join(" · ") || "조건 없음"];
    if (rule.ap) out.push(`원자 성질: ${apKo(rule.ap)}`);
    if (rule.ces && rule.ces.length) out.push(...cesKo(rule.ces));
    return out;
  }

  const PAIR = { cc: "cd", ce: "cf", cg: "ch", pc: "pd", pe: "pf", nc: "nd", ne: "nf" };

  function gaffExplain(r, label, idx) {
    const ff = r.ff, info = r.info, rules = ff.def.rules, g = r.graph;
    const c56 = new Set(info.c56);
    const groups = new Map();
    for (const i of idx) {
      const k = info.ruleIdx[i], rule = rules[k], pre = info.preTypes[i];
      const lines = gaffRuleLines(rule);
      const before = [...new Set(rules.slice(0, k).filter(x => x.z === rule.z && x.type !== pre).map(x => x.type))];
      if (before.length) lines.push(`이 규칙보다 앞에 있는 같은 원소 규칙(${before.join(", ")})에는 맞지 않음`);
      if (c56.has(i)) {
        lines.push(`규칙 결과는 ${pre} 이지만, moltemplate gaff2.lt 는 GAFF2 2.11 판이라 c5·c6 이 없어 2.11 규칙대로 c3`);
      } else if (pre in PAIR) {
        const a = pre, b = PAIR[pre];
        lines.push(`${a}/${b} 는 같은 원자의 두 이름: 공액 사슬에서 단일결합 양끝은 같은 이름(${a}–${a}), 이중·삼중결합 ` +
          `양끝은 다른 이름(${a}–${b})이 되도록 번갈아 붙임 — 결합 파라미터가 이 이름 쌍으로 정해짐 (antechamber atadjust)` +
          (pre !== label ? ` → 그래서 ${label}` : ` → ${label} 그대로`));
      } else if (pre === "cp") {
        lines.push("cp/cq 는 같은 원자의 두 이름: 연결 탄소끼리 단일결합이면 같은 이름, 아니면 다른 이름 (antechamber cpadjust)" +
          (pre !== label ? ` → 그래서 ${label}` : " → cp 그대로"));
      }
      const key = `${k}|${pre}|${c56.has(i)}`;
      if (!groups.has(key)) {
        groups.set(key, { title: `ATOMTYPE_GFF2.DEF ${k + 1}번째 규칙 — 위에서부터 처음 맞은 규칙`, code: rule.line || null, lines, atoms: [] });
      }
      groups.get(key).atoms.push(i);
    }
    const atoms = idx.map(i => {
      const lines = [];
      if (g.atoms[i].sym === "H") return { i, name: atomName(g, i), flagged: r.flagged.includes(i), env: environment(r, i), lines };
      const rg = info.rg[i], ar = info.ar[i];
      const sizes = [];
      for (let s = 3; s <= 10; s++) if (rg[s]) sizes.push(`${s}원${rg[s] > 1 ? ` ${rg[s]}개` : ""}`);
      const arc = [1, 2, 3, 4, 5].filter(k => ar[k] > 0).map(k => `AR${k}`);
      lines.push(`antechamber 고리 판정: ${sizes.length ? `${sizes.join(", ")} 고리${arc.length ? ` (${arc.join(", ")})` : ""}` : "고리 밖 (NR)"}`);
      const gb = info.gbonds.filter(([a, b]) => (a === i || b === i) && g.atoms[a === i ? b : a].sym !== "H")
        .map(([a, b, t]) => `${atomName(g, a === i ? b : a)} ${GBOND_KO[t] || t}`);
      if (gb.length) lines.push(`antechamber 결합 형식: ${gb.join(", ")}`);
      return { i, name: atomName(g, i), flagged: r.flagged.includes(i), env: environment(r, i), lines };
    });
    const mlines = [];
    if (GAFF2_KO[label]) mlines.push(GAFF2_KO[label]);
    mlines.push(`질량 ${ff.rules.masses[label]}`);
    return { meaning: { original: ff.desc[label] || null, source: "gaff2.lt", lines: mlines }, groups: [...groups.values()], atoms };
  }

  // ---------------------------------------------------------------- OPLS-AA
  function flatten(e, op) { return e[0] === op ? flatten(e[1], op).concat(flatten(e[2], op)) : [e]; }

  /** foyer SMARTS 원자 식 → 한국어. */
  function smartsExprKo(e) {
    switch (e[0]) {
      case "and": return flatten(e, "and").map(x => (x[0] === "or" ? `(${smartsExprKo(x)})` : smartsExprKo(x))).join(", ");
      case "or": return flatten(e, "or").map(x => (x[0] === "and" ? `(${smartsExprKo(x)})` : smartsExprKo(x))).join(" 또는 ");
      case "not": {
        const x = e[1];
        if (x[0] === "label") return `${x[1]} 규칙에 맞지 않은 원자`;
        if (x[0] === "r") return `${x[1]}원 고리 원자 아님`;
        if (x[0] === "X") return `이웃 ${x[1]}개 아님`;
        if (x[0] === "sym") return x[1] === "*" ? "원자 없음" : `${x[1]} 아님`;
        return `(${smartsExprKo(x)}) 아님`;
      }
      case "sym": return e[1] === "*" ? "아무 원자" : e[1];
      case "num": return `원자번호 ${e[1]} (${T.SYMBOLS[e[1]]})`;
      case "X": return `이웃 ${e[1]}개`;
      case "r": return `${e[1]}원 고리 원자`;
      case "R": return `고리 ${e[1]}개에 속함`;
      case "label": return `${e[1]} 규칙에 맞은 원자`;
      default: return JSON.stringify(e);
    }
  }

  /** 패턴을 나무 모양으로: 앞쪽 n-1 개 결합이 문자열의 사슬·가지, 나머지가 고리 닫기 (foyer_smarts 해석 순서). */
  function smartsTreeKo(pat) {
    const n = pat.atoms.length;
    const kids = Array.from({ length: n }, () => []);
    pat.edges.slice(0, n - 1).forEach(([p, c]) => kids[p].push(c));
    const lines = [];
    const walk = (p, d) => {
      lines.push(`${"　".repeat(d)}${d ? "└ " : ""}[${p + 1}] ${p === 0 ? "이 원자: " : ""}${smartsExprKo(pat.atoms[p])}`);
      for (const q of kids[p]) walk(q, d + 1);
    };
    walk(0, 0);
    for (const [a, b] of pat.edges.slice(n - 1)) lines.push(`고리 닫기: [${a + 1}] 과 [${b + 1}] 도 결합`);
    return lines;
  }

  function oplsExplain(r, label, idx) {
    const ff = r.ff, info = r.info, g = r.graph;
    const rulesByName = new Map(ff.foyerRules.map(x => [x.name, x]));
    const groups = new Map();
    for (const i of idx) {
      const fn = info.foyer[i], rule = rulesByName.get(fn), paired = info.mapped[i] !== label;
      const key = `${fn}|${paired}`;
      if (groups.has(key)) { groups.get(key).atoms.push(i); continue; }
      const lines = [];
      if (ff.foyerDesc[fn]) lines.push(`foyer 규칙 설명 (원문): ${ff.foyerDesc[fn]}`);
      lines.push("패턴 — [1] 이 이 타입을 받는 원자이고, 나머지는 그 둘레에 있어야 하는 원자");
      lines.push({ tree: smartsTreeKo(rule.pat) });
      if (rule.pat.atoms.length >= 3) lines.push("패턴에서 결합으로 이어지지 않은 두 원자는 분자에서도 결합하지 않아야 함 (foyer 맞춤 방식)");
      if (rule.overrides.length) lines.push(`이 규칙이 우선하는 타입 (함께 맞아도 이 타입이 이김): ${rule.overrides.join(", ")}`);
      lines.push(`foyer ${fn} → oplsaa2024.lt 의 같은 번호 @atom:${info.mapped[i]}`);
      if (paired) {
        const pr = ff.pairs.find(([a, , b]) => a === info.mapped[i] && b === label);
        lines.push(`OPLS/2020 짝: 이웃에 ${pr ? pr[1].join("·") : "?"} 타입이 있어 ${info.mapped[i]} 대신 ${label} ` +
          "(oplsaa2024.lt 가 OPLS/2020 전하로 바꾸며 둔 짝 — 그대로 두면 전하 합이 맞지 않음)");
      }
      groups.set(key, { title: `foyer 규칙 ${fn}`, code: ff.foyerDef[fn] || null, lines, atoms: [i] });
    }
    const atoms = idx.map(i => {
      const lines = [];
      const fn = info.foyer[i], white = [...info.white[i]].sort(), black = info.black[i];
      const others = white.filter(x => x !== fn).map(x => {
        const by = white.filter(y => (rulesByName.get(y) || { overrides: [] }).overrides.includes(x));
        return black.has(x) ? `${x} (${by.length ? by.join("·") + " 에 밀림" : "밀림"})` : x;
      });
      lines.push(`맞은 foyer 규칙: ${fn} (최종)` + (others.length ? `, ${others.join(", ")}` : ""));
      if (g.atoms[i].sym === "H") return { i, name: atomName(g, i), flagged: r.flagged.includes(i), env: environment(r, i), lines };
      const sizes = [...(info.cycles[i] || [])].map(c => c.split(",").length);
      lines.push(sizes.length ? `foyer 고리 판정: ${[...new Set(sizes)].sort((a, b) => a - b).map(s => `${s}원`).join(", ")} ` +
        `고리 (R = ${sizes.length})` : "foyer 고리 판정: 고리 아님 (8원 이하 고리 없음)");
      return { i, name: atomName(g, i), flagged: r.flagged.includes(i), env: environment(r, i), lines };
    });
    const mlines = [];
    const m = /_b([^_]*)_a([^_]*)_d([^_]*)_i(.*)$/.exec(ff.rules.full[label] || "");
    if (m) mlines.push(`결합 파라미터를 찾는 이름: 결합 ${m[1]} · 각 ${m[2]} · 이면각 ${m[3]} · improper ${m[4]}`);
    mlines.push(`타입 전하 ${ff.rules.charges[label]} · 질량 ${ff.rules.masses[label]}`);
    return { meaning: { original: ff.desc[label] || null, source: "oplsaa2024.lt", lines: mlines }, groups: [...groups.values()], atoms };
  }

  // ---------------------------------------------------------------- 진입점
  function explainType(r, label) {
    const idx = r.labels.map((l, i) => (l === label ? i : -1)).filter(i => i >= 0);
    let out;
    if (r.ff.NAME === "gaff2") out = gaffExplain(r, label, idx);
    else if (r.ff.NAME === "oplsaa") out = oplsExplain(r, label, idx);
    else out = dreidingExplain(r, label, idx);
    return Object.assign({ label, ffLabel: r.ff.LABEL, atomIdx: idx }, out);
  }

  /** 설명 전체를 글 줄로 (테스트·복사용). */
  function toText(x) {
    const flat = ls => ls.flatMap(l => (typeof l === "string" ? [l] : l.tree));
    const out = [`${x.label} (${x.ffLabel}) — 원자 ${x.atomIdx.map(i => i + 1).join(", ")}`];
    if (x.meaning.original) out.push(`원문 (${x.meaning.source}): ${x.meaning.original}`);
    out.push(...x.meaning.lines);
    for (const g of x.groups) {
      out.push(`${g.title} — 원자 ${g.atoms.map(i => i + 1).join(", ")}`);
      if (g.code) out.push(g.code);
      out.push(...flat(g.lines));
    }
    for (const a of x.atoms) out.push(`${a.name}${a.flagged ? " (경고 대상)" : ""}`, ...envText(a.env), ...a.lines);
    return out;
  }

  return { explainType, toText, environment, envText, neighborText, smallestRings, smartsExprKo, smartsTreeKo,
    dreidingMeaning, GAFF2_KO, BOND_SYM };
});
