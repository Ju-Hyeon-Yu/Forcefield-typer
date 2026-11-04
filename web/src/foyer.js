// fftyper/foyer_smarts.py·oplsaa.py 를 옮긴 foyer 규칙 OPLS-AA 판정 (패턴은 Python 이 해석한 JSON 을 그대로 씀)
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory(require("./typer.js"), require("./ltrules.js"));
  else root.FoyerTyper = factory(root.DreidingTyper, root.LtRules);
})(typeof self !== "undefined" ? self : this, function (T, L) {
  "use strict";
  // 판정 방식 원본: foyer 1.2.0 (MIT, © 2015 Vanderbilt University)

  const { TypingError, SYMBOLS } = T;

  // ---------------------------------------------------------------- 분자 그래프 (이웃은 결합이 나온 순서)
  class Graph {
    constructor(z, bonds) {
      this.n = z.length;
      this.z = z.slice();
      this.nbrs = z.map(() => []);
      this.adj = z.map(() => new Set());
      for (const [i, j] of bonds) {
        if (this.adj[i].has(j)) continue;
        this.nbrs[i].push(j); this.nbrs[j].push(i);
        this.adj[i].add(j); this.adj[j].add(i);
      }
      this.cycles = null;
    }

    /** foyer _find_chordless_cycles(max_cycle_size=8): 원자별 고리 튜플 집합 (문자열 키). */
    ringInfo() {
      if (this.cycles === null) this.cycles = chordlessCycles(this, 8).map(cs => new Set(cs.map(c => c.join(","))));
      return this.cycles;
    }
  }

  function chordlessCycles(g, maxSize) {
    const cycles = Array.from({ length: g.n }, () => []);
    for (let i = 0; i < g.n; i++) {
      const nb = g.nbrs[i];
      for (let p = 0; p < nb.length; p++) {
        for (let q = p + 1; q < nb.length; q++) {
          const a = nb[p], b = nb[q], last = a;
          let connected = false;
          let possible = [[a, i, b]];
          if (g.adj[last].has(b)) { cycles[i].push([a, i, b]); connected = true; }
          while (!connected) {
            const next = [];
            for (const ring of possible) {
              for (const x of g.nbrs[ring[ring.length - 1]]) if (x !== ring[ring.length - 2]) next.push(ring.concat([x]));
            }
            possible = next;
            for (const ring of possible) {
              const end = ring[ring.length - 1];
              if (g.adj[end].has(last)) {
                let chord = false;
                for (let k = 1; k < ring.length - 2; k++) if (g.adj[end].has(ring[k])) { chord = true; break; }
                if (!chord) { cycles[i].push(ring); connected = true; }
              }
            }
            if (!possible.length || possible[0].length === maxSize) break;
          }
        }
      }
    }
    return cycles;
  }

  // ---------------------------------------------------------------- 짝 맞추기 (유도 부분그래프)
  const ZOF = Object.fromEntries(SYMBOLS.map((s, z) => [s, z]).filter(([s]) => s));

  function atomOk(e, g, a, white) {
    switch (e[0]) {
      case "and": return atomOk(e[1], g, a, white) && atomOk(e[2], g, a, white);
      case "or": return atomOk(e[1], g, a, white) || atomOk(e[2], g, a, white);
      case "not": return !atomOk(e[1], g, a, white);
      case "sym": return e[1] === "*" || g.z[a] === ZOF[e[1]];
      case "num": return g.z[a] === e[1];
      case "label": return white[a].has(e[1]);
      case "X": return g.nbrs[a].length === e[1];
      case "r": for (const c of g.ringInfo()[a]) if (c.split(",").length === e[1]) return true; return false;
      case "R": return g.ringInfo()[a].size === e[1];
      default: throw new Error("알 수 없는 식 " + e[0]);
    }
  }

  function plan(pat) {
    const n = pat.atoms.length;
    const padj = Array.from({ length: n }, () => new Set());
    for (const [i, j] of pat.edges) { padj[i].add(j); padj[j].add(i); }
    const order = [[0, null]], seen = new Set([0]);
    for (let k = 0; k < order.length; k++) {
      const p = order[k][0];
      for (const q of [...padj[p]].sort((x, y) => x - y)) if (!seen.has(q)) { seen.add(q); order.push([q, p]); }
    }
    if (order.length !== n) throw new Error("이어지지 않은 패턴");
    return [order, padj];
  }

  function matchesAt(pat, g, root, white) {
    if (!pat._plan) pat._plan = plan(pat);
    const [order, padj] = pat._plan;
    const atoms = pat.atoms;
    if (!atomOk(atoms[0], g, root, white)) return false;
    const n = order.length;
    const image = new Map([[0, root]]);
    const used = new Set([root]);
    const rec = k => {
      if (k === n) return true;
      const [p, anchor] = order[k];
      for (const h of g.nbrs[image.get(anchor)]) {
        if (used.has(h) || !atomOk(atoms[p], g, h, white)) continue;
        let ok = true;
        for (const [q, hq] of image) if (padj[p].has(q) !== g.adj[h].has(hq)) { ok = false; break; }
        if (!ok) continue;
        image.set(p, h); used.add(h);
        if (rec(k + 1)) return true;
        image.delete(p); used.delete(h);
      }
      return false;
    };
    return rec(1);
  }

  /** foyer find_atomtypes 와 같은 순서 → 원자별 남은 후보 (정렬된 배열). */
  function findAtomtypes(g, rules, maxIter = 10) {
    const present = new Set(g.z.filter(z => z > 0).map(z => SYMBOLS[z]));
    const rs = rules.filter(r => r.pat.element === null || present.has(r.pat.element));
    const white = Array.from({ length: g.n }, () => new Set());
    const black = Array.from({ length: g.n }, () => new Set());
    for (let it = 0; it < maxIter; it++) {
      let found = false;
      for (const r of rs) {
        for (let a = 0; a < g.n; a++) {
          if (matchesAt(r.pat, g, a, white) && !white[a].has(r.name)) {
            white[a].add(r.name);
            for (const o of r.overrides) black[a].add(o);
            found = true;
          }
        }
      }
      if (!found) break;
    }
    return white.map((w, a) => [...w].filter(x => !black[a].has(x)).sort());
  }

  // ---------------------------------------------------------------- OPLS-AA (Python OplsTyper)
  class Opls {
    constructor(tables) {
      this.NAME = "oplsaa"; this.LABEL = "OPLS-AA"; this.LT_FILE = "oplsaa2024.lt"; this.LT_OBJECT = "OPLSAA";
      this.TYPE_CHARGES = true;
      this.source = tables.source;
      this.moltemplate = tables.moltemplate;
      this.foyerRules = tables.rules.map(r => ({ name: r.name, overrides: r.overrides, pat: r.pat }));
      this.foyerDesc = Object.fromEntries(tables.rules.map(r => [r.name, r.desc || ""]));
      this.no2024 = tables.no2024;
      this.pairs = tables.pairs;
      this.desc = tables.desc || {};
      this.rules = new L.LtRulesTable(tables.lt);
    }

    to2024(name) {
      const num = name.startsWith("opls_") ? name.slice(5) : "";
      if (name in this.no2024 || !/^\d+$/.test(num) || !(num in this.rules.charges)) return null;
      return num;
    }

    foyerTypes(g) {
      const fg = new Graph(g.atoms.map(a => SYMBOLS.indexOf(a.sym)), g.bonds.map(b => [b.a, b.b]));
      const cands = findAtomtypes(fg, this.foyerRules);
      const none = [], many = [];
      cands.forEach((c, i) => { if (!c.length) none.push(i); else if (c.length > 1) many.push(i); });
      if (none.length || many.length) {
        const sym = g.atoms.map(a => a.sym), msg = [];
        if (none.length) msg.push("OPLS-AA 규칙(foyer)에 맞는 타입이 없는 원자: " + none.map(i => `${i + 1}(${sym[i]})`).join(", "));
        if (many.length) msg.push("타입이 여럿 맞는 원자: " + many.map(i => `${i + 1}(${sym[i]}: ${cands[i].join(", ")})`).join(", "));
        throw new TypingError(msg.join("; "));
      }
      return cands.map(c => c[0]);
    }

    /** graph → {labels, warnings, flagged, extra} (Python OplsTyper.assign 과 같은 라벨·문장) */
    assign(g) {
      const names = this.foyerTypes(g);
      const labels = [], bad = new Map();
      names.forEach((nm, i) => {
        const num = this.to2024(nm);
        if (num === null || !(num in this.rules.masses)) {
          if (!bad.has(nm)) bad.set(nm, []);
          bad.get(nm).push(i + 1);
        }
        labels.push(num);
      });
      if (!bad.size) {
        const old = labels.slice();
        old.forEach((t, i) => {
          for (const [a, near, b] of this.pairs) {
            if (t === a && g.nbrs[i].some(([j]) => near.includes(old[j])) && b in this.rules.masses) labels[i] = b;
          }
        });
      }
      if (bad.size) {
        throw new TypingError("oplsaa2024.lt 에 대응하는 타입이 없습니다 — " + [...bad].map(([nm, at]) =>
          `원자 ${at.join(", ")}: foyer ${nm} (${this.foyerDesc[nm] || ""}) — ${this.no2024[nm] || "oplsaa2024.lt 에 같은 번호가 없음"}`)
          .join("; "));
      }
      const miss = L.missingTerms(g, labels, this.rules);
      return { labels, warnings: miss.warnings, flagged: miss.atoms, extra: "" };
    }

    charges(labels) { return labels.map(t => this.rules.charges[t]); }

    /** 화면용 설명: oplsaa2024.lt In Charges 주석의 설명 먼저, 그다음 결합 타입·전하. */
    describe(label) {
      const [left, text] = (this.desc[label] || "").split("|").map(x => (x || "").trim());
      const bt = (left || "").split("-").pop().trim();
      return `${text || left || "설명 없음"} (${bt ? bt + ", " : ""}전하 ${this.rules.charges[label]})`;
    }

    isMainLabel(label, sym) {
      if (sym !== "H") return true;
      const m = /_b([^_]*)_a/.exec(this.rules.full[label] || "");
      return !(m && (m[1] === "HC" || m[1] === "HA"));
    }
  }

  return { Graph, chordlessCycles, findAtomtypes, matchesAt, Opls };
});
