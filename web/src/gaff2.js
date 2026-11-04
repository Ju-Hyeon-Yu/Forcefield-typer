// fftyper/gaff2.py(AmberTools antechamber GAFF2 판정 이식)를 그대로 옮긴 브라우저/Node 용 GAFF2 타이퍼
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory(require("./typer.js"), require("./ltrules.js"));
  else root.Gaff2Typer = factory(root.DreidingTyper, root.LtRules);
})(typeof self !== "undefined" ? self : this, function (T, L) {
  "use strict";
  // 원본: AmberTools / AmberClassic antechamber (Junmei Wang 외), GPL-2.0-or-later.
  // 반복 순서·동률 처리를 Python 판과 같게 두어 결과가 원자 단위로 같아야 한다.

  const { TypingError, SYMBOLS } = T;

  // ---------------------------------------------------------------- 규칙 해석 (화학 환경의 원자 성질 문자열)
  function parseAp(text) {
    return text.replace(/^[\[\]]+|[\[\]]+$/g, "").split(",").map(cons => cons.split(".").map(unit => {
      const m = /^(\d*)([\s\S]*)$/.exec(unit);
      return [m[1] ? parseInt(m[1], 10) : -1, m[2]];
    }));
  }

  // ---------------------------------------------------------------- 분자 그래프 (antechamber ATOM/BOND 배열)
  class Mol {
    constructor(z, bonds) {
      this.n = z.length;
      this.z = z.slice();
      this.bonds = bonds.map(b => b.slice());
      this.con = z.map(() => []);
      for (const [i, j] of this.bonds) { this.con[i].push(j); this.con[j].push(i); }
      this.connum = this.con.map(c => c.length);
      this.byPair = new Map();
      for (const b of this.bonds) {
        const k1 = b[0] + "," + b[1], k2 = b[1] + "," + b[0];
        if (!this.byPair.has(k1)) this.byPair.set(k1, b);
        if (!this.byPair.has(k2)) this.byPair.set(k2, b);
      }
    }
    bondType(a, b) { const x = this.byPair.get(a + "," + b); return x ? x[2] : null; }
  }

  function ringOk(m, a) {
    const z = m.z[a], c = m.connum[a];
    if (![6, 7, 8, 15, 16].includes(z)) return false;
    if (z === 6 && c <= 2) return false;
    if ((z === 8 || z === 16) && c === 1) return false;
    return true;
  }

  function findRings(m) {
    const found = [], path = [];
    const cycle = start => {
      path.push(start);
      const sel = path.length;
      for (const nb of m.con[start].slice(0, 4)) {
        if (!ringOk(m, nb) || path.includes(nb)) continue;
        if (sel > 10) break;
        if (sel >= 2 && sel <= 9 && m.con[path[0]].slice(0, 4).includes(nb)) found.push(path.concat([nb]));
        cycle(nb);
      }
      path.pop();
    };
    for (let a = 0; a < m.n; a++) if (ringOk(m, a)) cycle(a);
    const rings = [], seen = new Set();
    for (const r of found) {
      const key = r.slice().sort((x, y) => x - y);
      const ks = key.join(",");
      if (seen.has(ks)) continue;
      seen.add(ks);
      const members = new Set(r);
      if (key.some(a => m.con[a].slice(0, 6).filter(x => members.has(x)).length === 3)) continue;
      rings.push(key);
    }
    return rings;
  }

  function ringdetect(m, rings, amberBcc) {
    const n = m.n;
    const rg = Array.from({ length: n }, () => new Array(11).fill(0));
    for (const r of rings) for (const a of r) { rg[a][0] += 1; rg[a][r.length] += 1; }
    const init = new Array(n).fill(0), ewd = new Array(n).fill(0);
    for (let a = 0; a < n; a++) {
      const z = m.z[a], c = m.connum[a];
      if (z === 6) init[a] = c === 3 ? 2 : (c === 4 ? -2 : 0);
      else if (z === 7) { ewd[a] = 1; init[a] = c <= 3 ? 2 : 0; }
      else if (z === 8) { ewd[a] = 1; init[a] = c === 2 ? 1 : 0; }
      else if (z === 15) init[a] = c === 2 ? 2 : (c === 3 ? 1 : (c >= 4 ? -1 : 0));
      else if (z === 16) { ewd[a] = 1; init[a] = c === 2 ? 1 : (c >= 3 ? -1 : 0); }
      else if ([9, 17, 35, 53].includes(z)) ewd[a] = 1;
    }
    const ar = Array.from({ length: n }, () => new Array(6).fill(0));
    for (const r of rings) {
      const size = r.length, tot = r.reduce((s, a) => s + init[a], 0);
      let kind = null;
      if (tot === -2 * size) kind = 5;
      else if (r.some(a => init[a] < 0)) kind = 4;
      else {
        if (size <= tot && tot <= 2 * size) {
          const inR = new Set(r);
          for (const [i, j, t] of m.bonds) {
            const idx = (inR.has(i) && rg[j][0] === 0 ? 1 : 0) + (inR.has(j) && rg[i][0] === 0 ? 1 : 0);
            if (idx === 1 && (t === 2 || t === 8)) { kind = 3; break; }
          }
        }
        if (kind === null && tot === 12 && size === 6) {
          let ok = true;
          for (const a of r) {
            if ((m.z[a] === 7 || m.z[a] === 15) &&
                !m.bonds.some(([i, j, t]) => (t === 8 || t === 2 || t === 10) && (i === a || j === a))) ok = false;
          }
          if (ok) kind = 1;
        }
        if (kind === null) kind = tot >= size + 3 ? 2 : 4;
      }
      for (const a of r) ar[a][kind] += 1;
    }
    if (amberBcc) {
      for (let a = 0; a < n; a++) {
        if (rg[a][5] >= 1 && rg[a][6] >= 1 && ar[a][1] + ar[a][2] >= 2) {
          for (const r of rings) if (r.length === 5 && r.includes(a)) for (const b of r) ar[b][2] -= 1;
        }
      }
    }
    const nr = ar.map(x => ([1, 2, 3, 4, 5].some(k => x[k] > 0) ? 0 : 1));
    return { rg, ar, nr, ewd };
  }

  const term = (m, a, z) => m.z[a] === z && m.connum[a] === 1;

  function bondtypePart(m, rings) {
    const { ar } = ringdetect(m, rings, true);
    const termOS = a => term(m, a, 8) || term(m, a, 16);
    const conj = new Array(m.n).fill(0);
    for (let a = 0; a < m.n; a++) {
      const z = m.z[a], c = m.connum[a], con = m.con[a];
      const cnt = k => con.slice(0, k).filter(termOS).length;
      if (z === 6 && c === 3 && cnt(3) >= 2) conj[a] = 1;
      if (z === 15 && c === 4 && cnt(4) >= 2) conj[a] = 1;
      if (z === 16 && c === 4 && cnt(4) >= 3) conj[a] = 1;
      if (z === 7 && c === 3 && cnt(3) >= 2) conj[a] = 1;
    }
    for (const b of m.bonds) {
      const bi = b[0], bj = b[1];
      if ((ar[bi][1] > 0 && ar[bj][1] > 0) || (ar[bi][2] > 0 && ar[bj][2] > 0)) {
        for (const r of rings) {
          if (r.length >= 7 || r.length <= 4) continue;
          if (r.some(a => ar[a][1] <= 0 && ar[a][2] <= 0)) continue;
          if (r.includes(bi) && r.includes(bj)) {
            if (b[2] === 1) b[2] = 7;
            else if (b[2] === 2) b[2] = 8;
            break;
          }
        }
      }
      if (b[2] === 2) {
        if (conj[bi] === 1 && termOS(bj)) { b[2] = 9; continue; }
        if (conj[bj] === 1 && termOS(bi)) { b[2] = 9; continue; }
      }
      if (m.z[bi] === 7 && (m.connum[bi] === 2 || m.connum[bi] === 3) && termOS(bj)) {
        for (const x of m.con[bi].slice(0, 6)) {
          if (x === bj) continue;
          if (termOS(x)) b[2] = 6;
          break;
        }
        if (b[2] === 6) continue;
      }
      if (m.z[bj] === 7 && (m.connum[bj] === 2 || m.connum[bj] === 3) && termOS(bi)) { b[2] = 6; continue; }
      if (b[2] === 1 && (termOS(bi) || termOS(bj))) { b[2] = 9; continue; }
    }
    for (let a = 0; a < m.n; a++) {
      let sel = null;
      if (m.z[a] === 16 && m.connum[a] === 3) {
        sel = m.con[a].slice(0, 3).filter(termOS);
        if (sel.length !== 2) sel = null;
      }
      if (sel) for (const b of m.bonds) if ((b[0] === a && sel.includes(b[1])) || (b[1] === a && sel.includes(b[0]))) b[2] = 9;
      if (m.z[a] === 15 && m.connum[a] === 4) {
        const s2 = m.con[a].slice(0, 4).filter(termOS);
        if (s2.length >= 2) for (const b of m.bonds) if ((b[0] === a && s2.includes(b[1])) || (b[1] === a && s2.includes(b[0]))) b[2] = 9;
      }
    }
  }

  // ---------------------------------------------------------------- 원자 타입 판정
  const JBOND = { 1: [1, 7, 9, 10], "-1": [1, 9], 2: [2, 8, 9, 10], "-2": [2, 9], 3: [3], "-3": [3], 4: [7, 8, 10], 9: [9] };
  const BOND_PROP = { SB: ["SB", -1], sb: ["sb", 1], DB: ["DB", -2], db: ["db", 2], TB: ["TB", -3], tb: ["tb", 3],
    DL: ["DL", 9], AB: ["AB", 4] };
  const COUNT_ADD = { 1: ["sb", "SB"], 2: ["db", "DB"], 3: ["tb", "TB"], 7: ["AB", "sb"], 8: ["AB", "db"],
    9: ["sb", "SB", "DL"], 10: ["AB"] };

  class Typer {
    constructor(m, rules) {
      this.m = m;
      this.rules = rules;
      this.wild = rules.wild;
      const rings = findRings(m);
      bondtypePart(m, rings);
      const p = ringdetect(m, rings, false);
      this.rg = p.rg; this.ar = p.ar; this.nr = p.nr; this.ewd = p.ewd;
      this.cnt = {};
      for (const k of ["sb", "SB", "db", "DB", "tb", "TB", "AB", "DL"]) this.cnt[k] = new Array(m.n).fill(0);
      for (const [i, j, t] of m.bonds) for (const k of COUNT_ADD[t] || []) { this.cnt[k][i] += 1; this.cnt[k][j] += 1; }
      this.apCache = new Map();
    }

    jbond(a, b, kind) { const t = this.m.bondType(a, b); return t !== null && JBOND[kind].includes(t); }

    apcheck(a, pre, ap) {
      for (const alts of ap) {
        let suc = false;
        for (const [count, name] of alts) if (this.prop(a, pre, count, name)) { suc = true; break; }
        if (!suc) return false;
      }
      return true;
    }

    prop(a, pre, count, name) {
      const cmp = v => count === v || (count === -1 && v > 0);
      if (name === "RG") return cmp(this.rg[a][0]);
      let m = /^RG(\d+)$/.exec(name);
      if (m) return parseInt(m[1], 10) <= 10 && cmp(this.rg[a][parseInt(m[1], 10)]);
      if (name === "NR") return cmp(this.nr[a]);
      m = /^AR([1-5])$/.exec(name);
      if (m) return cmp(this.ar[a][parseInt(m[1], 10)]);
      const key = name.slice(0, 2);
      if (Object.prototype.hasOwnProperty.call(BOND_PROP, key)) {
        let q = name.slice(2, 3) === "'" ? 1 : 0;
        if (q === 1 && name.slice(3, 4) === "'") q = 2;
        const [counter, kind] = BOND_PROP[key];
        if (!cmp(this.cnt[counter][a])) return false;
        if (q === 0) return true;
        if (pre === -1) return false;
        return this.jbond(a, pre, kind) === (q === 1);
      }
      return false;
    }

    nameMatch(a, name) {
      const z = this.m.z[a];
      const sym = z > 0 ? SYMBOLS[z] : "";
      if (sym.length === name.length && sym === name) return true;
      return Object.prototype.hasOwnProperty.call(this.wild, name) && this.wild[name].includes(z);
    }

    parsedAp(text) {
      if (!this.apCache.has(text)) this.apCache.set(text, parseAp(text));
      return this.apCache.get(text);
    }

    ces(center, chains) {
      const maxchain = Math.max(...chains.map(c => c.length));
      const found = chains.map(() => []);
      const order = [];
      const path = [center];
      const rec = () => {
        const length = path.length - 1;
        chains.forEach((ch, k) => {
          if (length !== ch.length) return;
          for (let j = 1; j <= length; j++) {
            const a = path[j], node = ch[j - 1];
            if (node.connum && this.m.connum[a] !== node.connum) return;
            if (node.name.startsWith("EW") && this.ewd[a] !== 1) return;
            if (!this.nameMatch(a, node.name)) return;
            if (node.ap && !this.apcheck(a, path[j - 1], this.parsedAp(node.ap))) return;
          }
          order.push([k, path.slice(1)]);
          found[k].push(order.length - 1);
        });
        for (const nb of this.m.con[path[path.length - 1]].slice(0, 6)) {
          if (path.includes(nb)) continue;
          if (path.length > maxchain) return;
          path.push(nb);
          rec();
          path.pop();
        }
      };
      rec();
      if (found.some(f => !f.length)) return false;
      const sel = new Array(chains.length).fill(-1);
      const check = () => {
        for (let i = 0; i < chains.length; i++) {
          for (let j = i + 1; j < chains.length; j++) {
            const t1 = sel[i], t2 = sel[j];
            if (t1 === t2) return false;
            const a1 = order[t1][1], a2 = order[t2][1];
            const mn = Math.min(a1.length, a2.length);
            let allSame = true;
            for (let k = 0; k < mn; k++) if (a1[k] !== a2[k]) { allSame = false; break; }
            if (allSame) return false;
            for (let k = 0; k < mn; k++) {
              const same = chains[i][k].cesname === chains[j][k].cesname;
              if ((a1[k] === a2[k]) !== same) return false;
            }
          }
        }
        return true;
      };
      const dc = level => {
        for (let idx = 0; idx < order.length; idx++) {
          if (order[idx][0] !== level) continue;
          sel[level] = idx;
          if (level === chains.length - 1) { if (check()) return true; }
          else if (dc(level + 1)) return true;
        }
        return false;
      };
      return dc(0);
    }

    match(a, rule) {
      const m = this.m;
      if (rule.res === "&") return true;
      if (rule.res !== "*") return false;
      for (const key of ["z", "connum", "nh", "new", "ap", "ces"]) {
        if (!(key in rule)) return true;
        const v = rule[key];
        if (key === "z" && v !== m.z[a]) return false;
        if (key === "connum" && v !== null && v !== m.connum[a]) return false;
        if (key === "nh" && v !== null && v !== m.con[a].slice(0, 6).filter(x => m.z[x] === 1).length) return false;
        if (key === "new" && v !== null) {
          if (!m.con[a].length || v !== m.con[m.con[a][0]].slice(0, 6).filter(x => this.ewd[x] === 1).length) return false;
        }
        if (key === "ap" && v !== null && !this.apcheck(a, -1, v)) return false;
        if (key === "ces" && !this.ces(a, v)) return false;
      }
      return true;
    }

    run() {
      const types = [];
      for (let a = 0; a < this.m.n; a++) {
        let t = null;
        for (const rule of this.rules.rules) if (this.match(a, rule)) { t = rule.type; break; }
        types.push(t);
      }
      this.atadjust(types);
      this.cpadjust(types);
      return types;
    }

    atadjust(types) {
      const pair = { cc: "cd", ce: "cf", cg: "ch", pc: "pd", pe: "pf", nc: "nd", ne: "nf" };
      const n = this.m.n, idx1 = new Array(n).fill(0), idx2 = new Array(n).fill(0);
      let num = 0, first = true;
      for (let a = 0; a < n; a++) {
        if (types[a] in pair) {
          idx2[a] = 1;
          if (first) { idx1[a] = 1; first = false; }
          num += 1;
        }
      }
      if (num === 0) return;
      for (let it = 0; it < num - 1; it++) {
        let flag = 0;
        for (const [i, j, t] of this.m.bonds) {
          if (idx2[i] + idx2[j] !== 2) continue;
          if (flag === 0 && idx1[i] === 0 && idx1[j] === 0) idx1[i] = 1;
          if (idx1[i] === 0 && idx1[j] !== 0) {
            flag = 1;
            if (t === 1 || t === 7) idx1[i] = idx1[j];
            if (t === 2 || t === 8 || t === 3) idx1[i] = -idx1[j];
          }
          if (idx1[j] === 0 && idx1[i] !== 0) {
            flag = 1;
            if (t === 1 || t === 7) idx1[j] = idx1[i];
            if (t === 2 || t === 8 || t === 3) idx1[j] = -idx1[i];
          }
        }
      }
      for (let a = 0; a < n; a++) if (idx1[a] === -1 && types[a] in pair) types[a] = pair[types[a]];
    }

    cpadjust(types) {
      const n = this.m.n, idx1 = new Array(n).fill(0), idx2 = new Array(n).fill(0);
      let num = 0, first = true;
      for (let a = 0; a < n; a++) {
        if (types[a] === "cp") {
          idx2[a] = 1;
          if (first) { idx1[a] = 1; first = false; }
          num += 1;
        }
      }
      if (num === 0) return;
      for (let it = 0; it < num - 1; it++) {
        for (const [i, j, t] of this.m.bonds) {
          if (idx2[i] + idx2[j] !== 2) continue;
          if (idx1[i] === 0 && idx1[j] !== 0) idx1[i] = t === 1 ? idx1[j] : -idx1[j];
          if (idx1[j] === 0 && idx1[i] !== 0) idx1[j] = t === 1 ? idx1[i] : -idx1[i];
        }
      }
      for (let a = 0; a < n; a++) if (idx1[a] === -1 && types[a] === "cp") types[a] = "cq";
    }

    consistencyWarnings(types) {
      const sign = { cc: 1, ce: 1, nc: 1, ne: 1, pc: 1, pe: 1, cg: 1, cd: -1, cf: -1, nd: -1, nf: -1, pd: -1, pf: -1, ch: -1 };
      const out = [];
      for (const [i, j, t] of this.m.bonds) {
        const si = sign[types[i]] || 0, sj = sign[types[j]] || 0;
        let bad = (si === sj && si !== 0 && !(t === 1 || t === 7)) || (si * sj === -1 && (t === 1 || t === 7));
        const pi = types[i], pj = types[j];
        if (((pi === "cp" && pj === "cp") || (pi === "cq" && pj === "cq")) && t !== 1) bad = true;
        if (((pi === "cp" && pj === "cq") || (pi === "cq" && pj === "cp")) && t === 1) bad = true;
        if (bad) out.push(`원자 ${i + 1}(${pi})–${j + 1}(${pj}): antechamber 도 이 공액 타입 조합이 틀릴 수 있다고 경고하는 경우입니다`);
      }
      return out;
    }
  }

  /** 웹 그래프(typer.js) → [원자번호 목록, [[i, j, 케쿨레 차수]]] — 결합 순서·방향 그대로. */
  function kekuleGraph(g) {
    const z = g.atoms.map(a => SYMBOLS.indexOf(a.sym));
    return [z, g.bonds.map(b => [b.a, b.b, b.kek])];
  }

  function assignTypes(z, bonds, rules) {
    const t = new Typer(new Mol(z, bonds), rules);
    const types = t.run();
    return [types, t.consistencyWarnings(types)];
  }

  // ---------------------------------------------------------------- 대응 타입 파라미터 (parmchk2 EQUA)
  const EQUA_LOOP = { bond: [0, 1], angle: [1, 0, 2], dihedral: [1, 2, 0, 3], improper: [0, 1, 2, 3] };

  function* equaCombos(equa, kind, types, places = null) {
    const order = EQUA_LOOP[kind];
    const lists = order.map(p => [types[p]].concat(places === null || places.includes(p) ? (equa[types[p]] || []) : []));
    const idx = new Array(lists.length).fill(0);
    const same = types.join("\u0000");
    while (true) {
      const t = new Array(types.length);
      order.forEach((p, k) => { t[p] = lists[k][idx[k]]; });
      if (t.join("\u0000") !== same) yield t;
      let k = lists.length - 1;
      while (k >= 0) { idx[k] += 1; if (idx[k] < lists[k].length) break; idx[k] = 0; k -= 1; }
      if (k < 0) return;
    }
  }

  function equaDecide(equa, rules, kind, types) {
    if (kind === "bond" || kind === "angle") {
      if (rules.match(kind, types) !== null) return ["ok", null];
      for (const t of equaCombos(equa, kind, types)) { const i = rules.match(kind, t); if (i !== null) return ["fill", i]; }
      return ["miss", null];
    }
    if (rules.match(kind, types, true, false) !== null) return ["ok", null];
    for (const t of equaCombos(equa, kind, types)) { const i = rules.match(kind, t, true, false); if (i !== null) return ["fill", i]; }
    if (rules.match(kind, types, false, true) !== null) return ["ok", null];
    for (const t of equaCombos(equa, kind, types, kind === "dihedral" ? [1, 2] : null)) {
      const i = rules.match(kind, t, false, true);
      if (i !== null) return ["fill", i];
    }
    return ["miss", null];
  }

  function canonical(kind, t, center = 2) {
    if (kind === "improper") {
      const outer = t.filter((_, i) => i !== center).sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
      outer.splice(center, 0, t[center]);
      return outer;
    }
    const rev = t.slice().reverse();
    if (kind === "dihedral") {
      return L.cmpTuple([t[1], t[2], t[0], t[3]], [t[2], t[1], t[3], t[0]]) <= 0 ? t.slice() : rev;
    }
    return L.cmpTuple(t, rev) <= 0 ? t.slice() : rev;
  }

  /** [[종류, 타입 배열, 가져온 계수 타입 이름]] — Python equa_fill 과 같은 순서. */
  function equaFill(g, labels, rules, equa) {
    const terms = L.topology(g, rules.improperCenter);
    const out = [], seen = new Set();
    for (const kind of L.KINDS) {
      if (kind === "improper" && rules.improperCenter === null) continue;
      for (const t of terms[kind]) {
        const types = t.map(i => labels[i]);
        const key = canonical(kind, types, rules.improperCenter);
        const sk = kind + "\u0001" + key.join("\u0000");
        if (seen.has(sk)) continue;
        seen.add(sk);
        const [what, i] = equaDecide(equa, rules, kind, types);
        if (what === "fill") {
          const src = rules.rules[kind][i][0];
          if (Object.prototype.hasOwnProperty.call(rules.coeffs[kind], src)) out.push([kind, key, src]);
        }
      }
    }
    return out;
  }

  function equaBlock(rules, adds, obj = "GAFF2") {
    if (!adds.length) return "";
    const names = new Map();
    const lines = ["# GAFF2 2.11 에 파라미터가 없는 항 — parmchk2 처럼 대응 타입(EQUA)의 gaff2.lt 값을 그대로 복사",
      `${obj} {`, '  write_once("In Settings") {'];
    for (const [kind, key, src] of adds) {
      let name = key.join("-");
      if (Object.prototype.hasOwnProperty.call(rules.coeffs[kind], name)) name += "_equa";
      names.set(kind + "\u0001" + key.join("\u0000"), name);
      lines.push(`    ${kind}_coeff @${kind}:${name} ${rules.coeffs[kind][src]}   # = ${src}`);
    }
    lines.push("  }");
    for (const kind of L.KINDS) {
      const rows = adds.filter(([k]) => k === kind);
      if (!rows.length) continue;
      lines.push(`  write_once("${rules.sections[kind]}") {`);
      for (const [, key] of rows) {
        lines.push(`    @${kind}:${names.get(kind + "\u0001" + key.join("\u0000"))} ` + key.map(x => `@atom:${x}`).join(" "));
      }
      lines.push("  }");
    }
    lines.push("}");
    return lines.join("\n") + "\n";
  }

  // ---------------------------------------------------------------- moltemplate gaff2.lt 용 타이퍼 (Python Gaff2Typer)
  class Gaff2 {
    constructor(tables) {
      this.NAME = "gaff2"; this.LABEL = "GAFF2"; this.LT_FILE = "gaff2.lt"; this.LT_OBJECT = "GAFF2";
      this.TYPE_CHARGES = false;
      this.def = tables.def;
      this.equa = tables.equa;
      this.desc = tables.desc || {};
      this.source = tables.source;
      this.moltemplate = tables.moltemplate;
      this.rules = new L.LtRulesTable(tables.lt);
    }

    /** graph → {labels, warnings, flagged, extra} */
    assign(g) {
      const [z, bonds] = kekuleGraph(g);
      const [types, warnings] = assignTypes(z, bonds, this.def);
      const masses = this.rules.masses;
      const old = [];
      types.forEach((t, i) => { if ((t === "c5" || t === "c6") && !(t in masses)) old.push(i); });
      for (const i of old) types[i] = "c3";
      if (old.length) {
        warnings.push(`원자 ${old.map(i => i + 1).join(", ")}: 이 gaff2.lt 는 GAFF2 2.11 판이라 ` +
          "5·6원 고리 sp3 C 를 c3 로 둡니다 (AmberTools GAFF2 2.2 의 c5·c6)");
      }
      const missing = [...new Set(types.filter(t => !(t in masses)))].sort();
      if (missing.length) {
        const where = [];
        types.forEach((t, i) => { if (missing.includes(t)) where.push(`${i + 1}(${t})`); });
        throw new TypingError(`gaff2.lt 에 없는 GAFF2 타입이 나왔습니다: ${missing.join(", ")} — 원자 ${where.join(", ")}`);
      }
      const adds = equaFill(g, types, this.rules, this.equa);
      if (adds.length) {
        warnings.push(`gaff2.lt 에 없는 항 ${adds.length}개는 대응 타입 값을 복사해 .lt 에 함께 씁니다 (parmchk2 EQUA): ` +
          adds.slice(0, 6).map(([k, key, src]) => `${L.KIND_KO[k]} ${key.join("-")}←${src}`).join(", ") +
          (adds.length > 6 ? " …" : ""));
      }
      const covered = {};
      for (const [kind, key] of adds) (covered[kind] = covered[kind] || new Set()).add(key.join("\u0000"));
      const miss = L.missingTerms(g, types, this.rules, covered);
      warnings.push(...miss.warnings);
      return { labels: types, warnings, flagged: miss.atoms, extra: equaBlock(this.rules, adds, this.LT_OBJECT) };
    }

    describe(label) {
      const m = this.rules.masses[label], d = this.desc[label];
      return `질량 ${m}${d ? " — " + d : ""}`;
    }

    isMainLabel(label, sym) { return sym !== "H" || ["hn", "ho", "hs", "hw"].includes(label); }
  }

  return { Gaff2, parseAp, findRings, kekuleGraph, assignTypes, equaFill, equaBlock, Typer, Mol };
});
