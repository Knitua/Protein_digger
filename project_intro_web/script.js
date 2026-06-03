const candidates = [
  {
    gene: "ALG13",
    accession: "Q9NP73",
    protein: "UDP-N-乙酰葡糖胺转移酶亚基 ALG13",
    deeploc: 0.5646,
    nls: 0.9678,
    ppi: 1,
    hiUnion: 1,
    rf2: 0,
    domain: "Tudor 结构域",
    strictDomain: null,
    qtm: 0.8798,
    lddt: 0.8441,
    pident: 16.0,
    matchedPlddt: 77.8255,
    annotation: "Tudor 结构域",
    annotationSourceShort: "UniProt / InterPro",
    annotationSource: "UniProt feature；InterPro IPR002999",
    interactionNote: "1 条调控锚点连接来自 HI-union 实验互作；RF2 final80 未提供额外高置信预测连接。",
    note: "已有 Tudor 结构域注释，更适合作为结构搜索流程回收已知/弱注释对象的阳性复核。"
  },
  {
    gene: "C11orf16",
    accession: "Q9NQ32",
    protein: "未表征蛋白 C11orf16",
    deeploc: 0.5662,
    nls: 0.9932,
    ppi: 1,
    hiUnion: 1,
    rf2: 0,
    domain: "PWWP 结构域；Tudor 结构域",
    strictDomain: "PWWP 结构域；Tudor 结构域",
    qtm: 0.8702,
    lddt: 0.8053,
    pident: 9.5,
    matchedPlddt: 75.7359,
    annotation: null,
    focus: true,
    interactionNote: "1 条调控锚点连接来自 HI-union 实验互作；提示该未表征蛋白与已知调控网络存在可追溯连接。",
    note: "DUF4537 区间被 PWWP 与 Tudor 双命中，是当前最值得进行人工结构叠合和功能假说复核的探索候选。"
  },
  {
    gene: "FBXO21",
    accession: "O94952",
    protein: "F-box 蛋白 21",
    deeploc: 0.5054,
    nls: 0.978,
    ppi: 1,
    hiUnion: 0,
    rf2: 1,
    domain: "Tudor 结构域",
    strictDomain: "Tudor 结构域",
    qtm: 0.8211,
    lddt: 0.7876,
    pident: 8.0,
    matchedPlddt: 91.1549,
    annotation: null,
    focus: true,
    interactionNote: "1 条调控锚点连接来自 RF2-PPI final80 高置信预测；当前缺少 HI-union 实验互作支持，后续需要优先复核预测界面和文献背景。",
    note: "Tudor 命中落在半甲基化 DNA 结合 / YccV-like 区间，结构可信但更像已知 DNA 结合折叠的交叉命中。"
  },
  {
    gene: "PHF21B",
    accession: "Q96EK2",
    protein: "PHD 指蛋白 21B",
    deeploc: 0.903,
    nls: 0.9995,
    ppi: 4,
    hiUnion: 3,
    rf2: 1,
    domain: "PHD 指结构域",
    strictDomain: null,
    qtm: 0.8789,
    lddt: 0.8431,
    pident: 41.1,
    matchedPlddt: 88.711,
    annotation: "PHD 指结构域",
    annotationSourceShort: "UniProt / Pfam / InterPro",
    annotationSource: "UniProt feature；Pfam PF00628；InterPro PHD zinc finger 条目",
    interactionNote: "共连接 4 个调控锚点，其中 3 条来自 HI-union 实验互作，1 条来自 RF2-PPI final80 高置信预测；互作证据强于多数候选。",
    note: "核定位与互作连接较强，但已有 PHD 指结构域注释，应作为调控相关候选复核，而非新结构域发现叙述。"
  },
  {
    gene: "PHF3",
    accession: "Q92576",
    protein: "PHD 指蛋白 3",
    deeploc: 0.8981,
    nls: 0.98,
    ppi: 1,
    hiUnion: 0,
    rf2: 1,
    domain: "PHD 指结构域",
    strictDomain: null,
    qtm: 0.8868,
    lddt: 0.8369,
    pident: 33.9,
    matchedPlddt: 88.6786,
    annotation: "PHD 指结构域",
    annotationSourceShort: "UniProt / Pfam / InterPro",
    annotationSource: "UniProt feature；Pfam PF00628；InterPro PHD zinc finger 条目",
    interactionNote: "1 条调控锚点连接来自 RF2-PPI final80 高置信预测；当前缺少 HI-union 实验互作支持。",
    note: "已有 PHD 指结构域注释，适合作为结构搜索流程回收已知同类结构域对象的验证项。"
  },
  {
    gene: "SMYD5",
    accession: "Q6GMV2",
    protein: "赖氨酸 N-三甲基转移酶 SMYD5",
    deeploc: 0.6522,
    nls: 0.9995,
    ppi: 2,
    hiUnion: 2,
    rf2: 0,
    domain: "SET 结构域",
    strictDomain: null,
    qtm: 0.8077,
    lddt: 0.7344,
    pident: 5.8,
    matchedPlddt: 93.4599,
    annotation: "SET 结构域",
    annotationSourceShort: "UniProt / Pfam / InterPro",
    annotationSource: "UniProt feature；Pfam PF00856；InterPro IPR001214 / IPR046341 / IPR044422",
    interactionNote: "共连接 2 个调控锚点，均来自 HI-union 实验互作；互作来源明确，但结构层面属于已知 SET 结构域回收。",
    note: "已有 SET 结构域注释，表观调控相关性明确，但不应被表述为未注释结构域发现。"
  },
  {
    gene: "ZCWPW2",
    accession: "Q504Y3",
    protein: "CW 型锌指 PWWP 结构域蛋白 2",
    deeploc: 0.9653,
    nls: 0.9512,
    ppi: 1,
    hiUnion: 1,
    rf2: 0,
    domain: "PWWP 结构域；Tudor 结构域",
    strictDomain: "Tudor 结构域",
    qtm: 0.8748,
    lddt: 0.8655,
    pident: 25.7,
    matchedPlddt: 93.0245,
    annotation: "PWWP 结构域",
    annotationSourceShort: "UniProt / Pfam / InterPro",
    focus: true,
    annotationSource: "UniProt feature；Pfam PF00855；InterPro IPR000313",
    interactionNote: "1 条调控锚点连接来自 HI-union 实验互作；结构上已有 PWWP 注释，但 Tudor 命中仍需要单独复核。",
    note: "已有 PWWP 注释，但 Tudor 命中未被同类注释过滤剔除，适合作为读码器样折叠交叉相似性的重点复核对象。"
  }
];

const listNode = document.querySelector("#candidateList");
const detailNode = document.querySelector("#candidateDetail");
const filterButtons = Array.from(document.querySelectorAll(".segment"));
let activeFilter = "all";
let selectedGene = "C11orf16";

function candidateMatches(candidate) {
  if (activeFilter === "strict") return Boolean(candidate.strictDomain);
  if (activeFilter === "review") return Boolean(candidate.focus);
  if (activeFilter === "annotated") return Boolean(candidate.annotation);
  return true;
}

function tagMarkup(candidate) {
  const tags = [`<span class="tag">${candidate.domain}</span>`];

  if (candidate.focus) tags.push(`<span class="tag review">重点人工复核</span>`);
  if (candidate.annotation) {
    tags.push(`<span class="tag warn">已有同类注释：${candidate.annotation}（${candidate.annotationSourceShort}）</span>`);
  } else {
    tags.push(`<span class="tag">未见同类注释</span>`);
  }
  return tags.join("");
}

function renderCandidates() {
  const filtered = candidates.filter(candidateMatches);
  if (!filtered.some((candidate) => candidate.gene === selectedGene)) {
    selectedGene = filtered[0]?.gene || candidates[0].gene;
  }

  listNode.innerHTML = filtered
    .map(
      (candidate) => `
        <article class="candidate-card ${candidate.gene === selectedGene ? "active" : ""}" data-gene="${candidate.gene}" tabindex="0">
          <div class="candidate-title">
            <h3>${candidate.gene}</h3>
            <span class="accession">${candidate.accession}</span>
          </div>
          <p class="candidate-protein">${candidate.protein}</p>
          <div class="tag-row">${tagMarkup(candidate)}</div>
        </article>
      `
    )
    .join("");

  listNode.querySelectorAll(".candidate-card").forEach((card) => {
    const select = () => {
      selectedGene = card.dataset.gene;
      renderCandidates();
      renderDetail();
    };
    card.addEventListener("click", select);
    card.addEventListener("keydown", (event) => {
      if (event.key === "Enter" || event.key === " ") {
        event.preventDefault();
        select();
      }
    });
  });
}

function scoreRow(label, value, max = 1) {
  const width = Math.max(0, Math.min(100, (value / max) * 100));
  return `
    <div class="score-row">
      <span>${label}</span>
      <span class="score-track"><span class="score-fill" style="width: ${width}%"></span></span>
      <strong>${value.toFixed(4)}</strong>
    </div>
  `;
}

function renderDetail() {
  const candidate = candidates.find((item) => item.gene === selectedGene) || candidates[0];
  const annotationText = candidate.annotation
    ? `${candidate.annotation}；来源：${candidate.annotationSource}`
    : "未见 UniProt feature / Pfam / InterPro 同类结构域注释";

  detailNode.innerHTML = `
    <h3>${candidate.gene}</h3>
    <p class="detail-protein">${candidate.protein}</p>
    <div class="tag-row">${tagMarkup(candidate)}</div>
    <div class="score-list">
      ${scoreRow("DeepLoc 核定位", candidate.deeploc)}
      ${scoreRow("NLSExplorer 最高分", candidate.nls)}
      ${scoreRow("Foldseek 最佳 qTM", candidate.qtm)}
      ${scoreRow("Foldseek 最佳 lDDT", candidate.lddt)}
      ${scoreRow("Foldseek pident（序列相似度 %）", candidate.pident, 100)}
      ${scoreRow("命中区间 mean pLDDT", candidate.matchedPlddt, 100)}
    </div>
    <div class="detail-note">
      <div class="detail-section">
        <strong>互作证据</strong>
        <span>与 ${candidate.ppi} 个 TF / EpiFactors 锚点相连；其中 HI-union 实验互作 ${candidate.hiUnion} 条，RF2-PPI final80 高置信预测 ${candidate.rf2} 条。${candidate.interactionNote}</span>
      </div>
      <div class="detail-section">
        <strong>结构证据</strong>
        <span>Foldseek 命中 ${candidate.domain}；最佳 qTM ${candidate.qtm.toFixed(4)}，最佳 lDDT ${candidate.lddt.toFixed(4)}，pident ${candidate.pident.toFixed(1)}%，target matched-region mean pLDDT ${candidate.matchedPlddt.toFixed(4)}。</span>
      </div>
      <div class="detail-section">
        <strong>注释来源</strong>
        <span>${annotationText}</span>
      </div>
      <div class="detail-section">
        <strong>复核解释</strong>
        <span>${candidate.note}</span>
      </div>
    </div>
  `;
}

filterButtons.forEach((button) => {
  button.addEventListener("click", () => {
    activeFilter = button.dataset.filter;
    filterButtons.forEach((item) => item.classList.toggle("active", item === button));
    renderCandidates();
    renderDetail();
  });
});

renderCandidates();
renderDetail();
