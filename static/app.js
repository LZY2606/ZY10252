const $ = (sel, root = document) => root.querySelector(sel);
const esc = (s) => String(s).replace(/[&<>"']/g, (c) =>
  ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

async function refresh() {
  const r = await fetch('/api/snapshot');
  const data = await r.json();
  render(data);
}

function postForm(url, params) {
  const body = new URLSearchParams(params);
  return fetch(url, { method: 'POST', body }).then(async (r) => {
    const j = await r.json();
    if (!r.ok || j.ok === false) throw new Error(j.error || '请求失败');
    return j;
  });
}

function atomLabel(a) {
  return a.isotope ? `${a.label.replace(/\[\d+\]$/, '')}${a.isotope ? '·' + a.isotope : ''}` : a.label;
}

// 极小的分子图：原子按圆周布局，键为直线。
function moleculeSvg(mol, size = 180) {
  const n = Math.max(mol.atoms.length, 1);
  const cx = size / 2, cy = size / 2, R = size * 0.34;
  const pos = mol.atoms.map((_, i) => {
    const ang = -Math.PI / 2 + (2 * Math.PI * i) / n;
    return [cx + R * Math.cos(ang), cy + R * Math.sin(ang)];
  });
  let edges = mol.bonds.map((b) => {
    const [x1, y1] = pos[b.u], [x2, y2] = pos[b.v];
    const off = b.order === 2 ? 2 : b.order === 3 ? 4 : 0;
    if (!off) return `<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke-width="1.6"/>`;
    const dx = x2 - x1, dy = y2 - y1, L = Math.hypot(dx, dy) || 1, nx = -dy / L * off, ny = dx / L * off;
    let out = `<line x1="${x1 + nx}" y1="${y1 + ny}" x2="${x2 + nx}" y2="${y2 + ny}" stroke-width="1.4"/>
               <line x1="${x1 - nx}" y1="${y1 - ny}" x2="${x2 - nx}" y2="${y2 - ny}" stroke-width="1.4"/>`;
    if (b.order === 3) out += `<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke-width="1.2"/>`;
    return out;
  }).join('');
  let nodes = mol.atoms.map((a, i) => {
    const [x, y] = pos[i];
    const base = a.label.replace(/\[\d+\]$/, '');
    const iso = a.isotope ? `<tspan font-size="8" baseline-shift="super">${a.isotope}</tspan>` : '';
    const ch = a.charge ? `<tspan font-size="9">${a.charge > 0 ? '+' : '−'}${Math.abs(a.charge)}</tspan>` : '';
    return `<circle cx="${x}" cy="${y}" r="13" fill="#0f1420" stroke="#60a5fa"/>
            <text x="${x}" y="${y + 3.6}" text-anchor="middle">${iso}${esc(base)}${ch}<tspan font-size="8" fill="#94a3b8">[${i}]</tspan></text>`;
  }).join('');
  return `<svg width="${size}" height="${size}" viewBox="0 0 ${size} ${size}">${edges}${nodes}</svg>`;
}

function render(data) {
  $('#fixture-meta').textContent = `fixture ${data.fixture} · sha256 ${data.sha256.slice(0, 16)}…`;
  const s = data.snapshot;

  $('#molecules').innerHTML = s.molecules.map((m) => `
    <div class="card">
      <h3>${esc(m.id)} <span class="tag ${m.charge === 0 ? 'ok' : 'warn'}">电荷 ${m.charge}</span></h3>
      ${m.atoms.length ? moleculeSvg(m) : '<div class="muted">抽象池物种（无组成，仅用于网络流量）</div>'}
    </div>`).join('');

  $('#reactions').innerHTML = s.reactions.map((r) => reactionCard(r)).join('');
  $('#networks').innerHTML = s.networks.map((n) => networkCard(n)).join('');
  bindReactionActions();

  $('#log-table tbody').innerHTML = s.log.map((e) => `
    <tr><td>${e.seq}</td><td>${esc(e.at)}</td><td>${esc(e.summary)}</td>
    <td class="muted edge">${esc(e.op_json)}</td></tr>`).join('');
}

function termList(terms) {
  return terms.map((t) =>
    `${t.coef === '1' ? '' : esc(t.coef) + ' '}<b>${esc(t.mol)}</b>${t.builtin ? '' : ' <span class="tag">自定义</span>'}`
  ).join(' &nbsp;+&nbsp; ');
}

function reactionCard(r) {
  const consRows = r.conservation.map((c) => `
    <tr><td>${esc(c.key)}</td>
      <td class="${c.balanced ? 'comp-bal' : 'comp-unbal'}">${esc(c.delta)}</td>
      <td>${c.balanced ? '✓ 守恒' : '✗ 未守恒'}</td></tr>`).join('');

  const mapClasses = r.map_classes.length ? r.map_classes.map((c) => {
    const locked = r.locked_class === c.class_id;
    const edges = c.edges.map((e) => `<div class="edge">${esc(e.from)} → ${esc(e.to)}</div>`).join('');
    return `<div class="mapclass ${locked ? 'locked' : ''}">
      <div class="row" style="justify-content:space-between">
        <div><b>等价类</b> <code>${c.class_id.slice(0, 12)}</code>
          <span class="tag ${c.multiplicity > 1 ? 'warn' : ''}">等价映射 ×${c.multiplicity}</span>
          ${locked ? '<span class="tag ok">已锁定</span>' : ''}
          <span class="tag ${c.identity_ok ? 'ok' : 'bad'}">身份${c.identity_ok ? '守恒' : '失败'}</span>
        </div>
        <div>
          ${locked
            ? `<button data-act="unlock" data-rxn="${esc(r.id)}">解锁</button>`
            : `<button class="primary" data-act="lock" data-rxn="${esc(r.id)}" data-class="${esc(c.class_id)}">锁定这组映射</button>`}
        </div>
      </div>
      ${edges}
      ${c.issues.length ? `<div class="comp-unbal">${c.issues.map(esc).join('；')}</div>` : ''}
    </div>`;
  }).join('') : '<div class="muted">该反应没有原子映射。</div>';

  const candidates = r.candidates.length ? r.candidates.map((c) => `
    <div class="cand">
      <div><b>候选修复</b> <span class="tag">证据代价 ${c.cost}</span></div>
      ${c.actions.map((a) => `<div class="edge">向${a.side === 'reactant' ? '反应物' : '产物'}补 ${esc(a.coef)} ${esc(a.mol)}</div>`).join('')}
      <div class="muted">${c.evidence.map(esc).join('；')}</div>
      <div class="section-mini">
        ${c.actions.map((a) => `<button data-act="add" data-rxn="${esc(r.id)}" data-side="${a.side}" data-mol="${esc(a.mol)}" data-coef="${esc(a.coef)}">应用：+${esc(a.coef)} ${esc(a.mol)}</button>`).join(' ')}
      </div>
    </div>`).join('')
    : (r.element_balanced && r.charge_balanced ? '<div class="muted">已守恒，无需修复。</div>' : '<div class="muted">有界枚举内未找到只靠补入内置小分子的修复。</div>');

  return `<div class="card section-mini">
    <h3>
      <span>${esc(r.id)}</span>
      <span>
        <span class="tag ${r.element_balanced ? 'ok' : 'bad'}">元素/同位素 ${r.element_balanced ? '守恒' : '未守恒'}</span>
        <span class="tag ${r.charge_balanced ? 'ok' : 'bad'}">电荷 ${r.charge_balanced ? '守恒' : '未守恒'}</span>
        <span class="tag ${r.direction_unknown ? 'warn' : ''}">${r.direction_unknown ? '方向未定' : '方向已定'}</span>
      </span>
    </h3>
    <div class="eqline">${termList(r.reactants)} &nbsp;⟶&nbsp; ${termList(r.products)}</div>
    <div class="row">
      <div class="col">
        <h3 style="font-size:13px">未守恒向量（反应物 − 产物，有理数）</h3>
        <table><thead><tr><th>分量</th><th>Δ</th><th>结论</th></tr></thead><tbody>${consRows}</tbody></table>
        <div class="section-mini">
          <button data-act="toggle-dir" data-rxn="${esc(r.id)}" data-unknown="${r.direction_unknown ? 'false' : 'true'}">
            ${r.direction_unknown ? '标记方向为已定' : '标记方向为未定'}
          </button>
        </div>
      </div>
      <div class="col">
        <h3 style="font-size:13px">原子映射边与对称等价类</h3>
        ${mapClasses}
        <h3 style="font-size:13px">有界候选修复（保留操作与证据代价）</h3>
        ${candidates}
      </div>
    </div>
  </div>`;
}

function networkCard(n) {
  const sols = n.solutions.length ? n.solutions.map((sol) => {
    const fl = sol.fluxes.map(([k, v]) => `f(${k}) = ${v}`).join(',  ');
    const inter = sol.intermediates.map(([k, v]) => `${k}: 净 ${v}`).join(',  ');
    return `<div class="cand"><div><span class="tag">分支</span> ${esc(sol.branch)}</div>
      <div class="sols">${esc(fl)}\n中间体稳态：${esc(inter)}</div></div>`;
  }).join('') : '<div class="muted">有界网格内无可行稳态流量。</div>';
  return `<div class="card">
    <h3>${esc(n.id)} <span class="tag">${n.steps.length} 步</span></h3>
    <div class="muted">步骤：${n.steps.map(esc).join(' → ')}</div>
    <div class="section-mini">${sols}</div>
  </div>`;
}

function bindReactionActions() {
  document.querySelectorAll('button[data-act]').forEach((btn) => {
    btn.onclick = async () => {
      try {
        const d = btn.dataset;
        if (d.act === 'lock') await postForm('/api/ops/lock', { reaction: d.rxn, class_id: d.class });
        if (d.act === 'unlock') await postForm('/api/ops/unlock', { reaction: d.rxn });
        if (d.act === 'add') await postForm('/api/ops/add-participant', { reaction: d.rxn, side: d.side, mol: d.mol, coef: d.coef });
        if (d.act === 'toggle-dir') await postForm('/api/ops/direction', { reaction: d.rxn, unknown: d.unknown });
        await refresh();
      } catch (e) { alert(e.message); }
    };
  });
}

$('#btn-refresh').onclick = refresh;
$('#btn-reset').onclick = async () => {
  if (!confirm('确定清空数据库并重新导入固定 fixture？')) return;
  await fetch('/api/reset', { method: 'POST' });
  await refresh();
};
$('#btn-export').onclick = async () => {
  const r = await fetch('/api/export');
  const blob = new Blob([JSON.stringify(await r.json(), null, 2)], { type: 'application/json' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = 'conservation-run.json';
  a.click();
};
$('#btn-reimport').onclick = async () => {
  const content = $('#reimport-text').value;
  const msg = $('#reimport-msg');
  try {
    const r = await fetch('/api/reimport', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ content, name: 'manual.rxn' }),
    });
    const j = await r.json();
    if (!r.ok) throw new Error(j.error || '导入失败');
    msg.textContent = '导入成功，已重放。';
    await refresh();
  } catch (e) { msg.textContent = '失败：' + e.message; }
};

refresh();
