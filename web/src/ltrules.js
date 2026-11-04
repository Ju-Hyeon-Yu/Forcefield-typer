// fftyper/ltrules.py 를 옮긴 moltemplate 역장 규칙표 — 결합·각·이면각·improper 규칙 찾기와 누락 경고
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.LtRules = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  const KINDS = ["bond", "angle", "dihedral", "improper"];
  const KIND_KO = { bond: "결합", angle: "결합각", dihedral: "이면각", improper: "improper" };

  /** Python 튜플 비교 (문자열 원소, 사전식). */
  function cmpTuple(a, b) {
    for (let k = 0; k < Math.min(a.length, b.length); k++) {
      if (a[k] < b[k]) return -1;
      if (a[k] > b[k]) return 1;
    }
    return a.length - b.length;
  }

  /** fnmatch.fnmatchcase 의 *, ?, [...] 만 지원하는 정규식. */
  function globRe(p) {
    let s = "";
    for (let i = 0; i < p.length; i++) {
      const c = p[i];
      if (c === "*") s += ".*";
      else if (c === "?") s += ".";
      else if (c === "[") {
        const j = p.indexOf("]", i + 1);
        if (j < 0) { s += "\\["; continue; }
        let body = p.slice(i + 1, j);
        if (body[0] === "!") body = "^" + body.slice(1);
        s += "[" + body.replace(/\\/g, "\\\\") + "]";
        i = j;
      } else s += c.replace(/[.+^${}()|\\]/g, "\\$&");
    }
    return new RegExp("^(?:" + s + ")$", "s");
  }

  function permutations(xs) {
    if (xs.length <= 1) return [xs.slice()];
    const out = [];
    xs.forEach((x, i) => {
      for (const p of permutations(xs.slice(0, i).concat(xs.slice(i + 1)))) out.push([x].concat(p));
    });
    return out;
  }

  function combinations(xs, k) {
    const out = [];
    const rec = (start, cur) => {
      if (cur.length === k) { out.push(cur.slice()); return; }
      for (let i = start; i < xs.length; i++) { cur.push(xs[i]); rec(i + 1, cur); cur.pop(); }
    };
    rec(0, []);
    return out;
  }

  /**
   * table: {masses, rules: {종류: [[이름, [패턴…]], …]}, coeffs: {종류: {이름: 계수}}, sections, improperCenter,
   *         full: {짧은 이름: replace{} 긴 이름}, charges: {타입: In Charges 전하}}
   * (web/build_tables.py 가 Python LtRules 에서 그대로 만든다.)
   */
  class LtRulesTable {
    constructor(table) {
      this.masses = table.masses;
      this.rules = table.rules;
      this.coeffs = table.coeffs;
      this.sections = table.sections;
      this.improperCenter = table.improperCenter === undefined ? null : table.improperCenter;
      this.full = table.full || {};
      this.charges = table.charges || {};
      this._exact = {};
      this._wild = {};
      for (const kind of KINDS) {
        const ex = new Map(), wild = [];
        (this.rules[kind] || []).forEach(([, pats], i) => {
          if (pats.some(p => /[*?[]/.test(p))) wild.push([i, pats.map(globRe)]);
          else ex.set(pats.join("\u0000"), i);
        });
        wild.reverse();
        this._exact[kind] = ex;
        this._wild[kind] = wild;
      }
    }

    orders(kind, names) {
      if (kind !== "improper") return [names.slice(), names.slice().reverse()];
      const c = this.improperCenter;
      const outer = names.filter((_, i) => i !== c);
      return permutations(outer).map(p => { const t = p.slice(); t.splice(c, 0, names[c]); return t; });
    }

    /** moltemplate 가 줄 규칙 번호(맞는 것 중 마지막) 또는 null. */
    match(kind, names, exact = true, wild = true) {
      let best = -1;
      names = names.map(x => (Object.prototype.hasOwnProperty.call(this.full, x) ? this.full[x] : x));
      for (const t of this.orders(kind, names)) {
        if (exact) {
          const i = this._exact[kind].get(t.join("\u0000"));
          if (i !== undefined && i > best) best = i;
        }
        if (wild) {
          for (const [i, res] of this._wild[kind]) {
            if (i <= best) break;
            if (res.every((re, k) => re.test(t[k]))) { best = i; break; }
          }
        }
      }
      return best >= 0 ? best : null;
    }
  }

  /** 그래프 → {종류: 원자 번호 배열 목록} (Python topology 와 같은 순서). */
  function topology(g, improperCenter = null) {
    const nb = g.nbrs.map(x => x.map(([j]) => j).sort((a, b) => a - b));
    const bonds = g.bonds.map(b => [Math.min(b.a, b.b), Math.max(b.a, b.b)])
      .sort((u, v) => u[0] - v[0] || u[1] - v[1]);
    const angles = [];
    nb.forEach((list, j) => list.forEach((i, x) => list.slice(x + 1).forEach(k => angles.push([i, j, k]))));
    const dihedrals = [];
    for (const [j, k] of bonds) {
      for (const i of nb[j]) {
        if (i === k) continue;
        for (const l of nb[k]) if (l !== j && l !== i) dihedrals.push([i, j, k, l]);
      }
    }
    const impropers = [];
    if (improperCenter !== null) {
      nb.forEach((list, c) => {
        for (const trio of combinations(list, 3)) { const t = trio.slice(); t.splice(improperCenter, 0, c); impropers.push(t); }
      });
    }
    return { bond: bonds, angle: angles, dihedral: dihedrals, improper: impropers };
  }

  /**
   * 파라미터 규칙이 없는 결합·각·이면각 경고 (Python missing_terms 와 같은 문장·순서).
   * covered: {종류: Set("a\u0000b…")} — .lt 에 따로 채워 넣는 항. 반환 {warnings, atoms(관련 원자)}.
   */
  function missingTerms(g, labels, rules, covered = {}, limit = 6) {
    const out = [], atoms = new Set();
    const terms = topology(g);
    for (const kind of ["bond", "angle", "dihedral"]) {
      const miss = new Map();
      const extra = covered[kind] || new Set();
      for (const t of terms[kind]) {
        const key = t.map(i => labels[i]), rev = key.slice().reverse();
        if (extra.has(key.join("\u0000")) || extra.has(rev.join("\u0000")) || rules.match(kind, key) !== null) continue;
        const k = cmpTuple(key, rev) <= 0 ? key : rev;
        const s = k.join("\u0000");
        if (!miss.has(s)) miss.set(s, [k, []]);
        miss.get(s)[1].push(t);
        t.forEach(i => atoms.add(i));
      }
      const items = [...miss.values()].sort((a, b) => cmpTuple(a[0], b[0]));
      for (const [key, where] of items) {
        const at = where.slice(0, 3).map(t => t.map(i => i + 1).join("-")).join(", ");
        const more = where.length > 3 ? ` 외 ${where.length - 3}개` : "";
        const what = kind === "bond" ? "moltemplate 가 오류로 멈춥니다 (No bond types defined)"
          : "moltemplate 가 이 항을 빼고 만듭니다";
        out.push(`${KIND_KO[kind]} ${key.join("-")} (원자 ${at}${more}): 역장 파일에 파라미터가 없어 ${what}`);
        if (out.length >= limit * 3) return { warnings: out, atoms: [...atoms].sort((a, b) => a - b) };
      }
    }
    return { warnings: out, atoms: [...atoms].sort((a, b) => a - b) };
  }

  return { KINDS, KIND_KO, LtRulesTable, topology, missingTerms, cmpTuple, globRe };
});
