export type Route = "A" | "B1" | "B2" | "B3";
export type Candidate = {id:string;row:number;route:Route;gene:string;accessions:string[];name:string};
export type Pair = {id:string;row:number;route:Route;gene:string;candidate:string;anchorGene:string;anchor:string;evidence:string;plm:number|null;pplm:number|null};
export type Release = {
 release:string;methodReview:string;
 counts:{excel_sha256:string;candidate_rows:number;unique_genes:number;unique_accessions:number;pair_rows:number;unique_pairs:number;appendix_rows:number;route_rows:Record<Route,number>;route_pair_rows:Record<Route,number>};
 candidates:Candidate[];pairs:Pair[];appendix:Record<string,string>[];
};
export type Doc = {id:string;file:string;title:string;section:string;sha256:string;characters:number};
export const routeRules:Record<Route,string> = {A:"功能证据分层 · HI-union ∪ RF2-final90",B1:"DeepLoc × NLSExplorer · HI-union ∪ RF2-final80",B2:"网络邻近度前 2% · HI-union ∪ RF2-final80",B3:"PLM-interact > 0.99 且 PPLM-PPI > 0.9"};
export const routeLabels:Record<Route,string> = {A:"核定位注释与功能分层",B1:"序列驱动的核定位预测",B2:"调控网络邻近性筛选",B3:"可变异构体核定位筛选"};
export const sections = [
 {id:"stage0",label:"蛋白组分区",title:"蛋白组分区与锚点排除",kicker:"代表序列的核定位注释分区",summary:"20,416 条人类审校代表序列，按适用于该序列的核定位注释分为 5,669 条与 14,747 条；排除调控锚点后，分别得到 3,418 条与 14,625 条。",version:"采用 2,368 个调控锚点，联合蛋白编号与独立基因符号排除。核定位注释必须适用于代表序列；可变异构体单独评价。"},
 {id:"a",label:"核定位注释与功能分层",title:routeLabels.A,kicker:"核定位证据 · 功能新颖性评价",summary:"3,418 条非锚点核定位蛋白，经四类调控证据与功能分层，保留 2,457 条聚焦候选；互作证据筛选后得到 954 条结果。",version:"四类一级证据依次为 541 / 236 / 440 / 2,201 条。聚焦集合由 440 条已知关联、924 条研究不足与 1,093 条混合功能候选组成。"},
 {id:"b1",label:"序列驱动的核定位预测",title:routeLabels.B1,kicker:"DeepLoc × NLSExplorer",summary:"14,625 条非核注释代表序列，经双模型核定位判定得到 1,350 条候选；HI-union 与 RF2-final80 并集支持其中 429 条。",version:"NLSExplorer 按相同序列长度组批，避免补齐位点影响推理。正式候选来自完整输入的一致判定，历史比较集合仅供核验。"},
 {id:"b2",label:"调控网络邻近性筛选",title:routeLabels.B2,kicker:"多源网络 · 调控种子邻近度",summary:"从 14,625 条输入中排除 1,350 条 B1 候选，对剩余 13,275 条排序；取前 2%（向上取整）得到 266 条，互作证据支持其中 85 条。",version:"使用 2,368 个调控锚点构图；四项指标先按最小名次处理并列并转换为百分位，再等权聚合。"},
 {id:"b3",label:"可变异构体核定位筛选",title:routeLabels.B3,kicker:"异构体特异性 · 联合互作预测",summary:"保守分析分支的 497 条异构体，与 2,368 个调控锚点组合评价；490 条配对通过联合条件，涉及 152 条异构体、112 个基因。",version:"本批保守分支使用历史 B1 的 1,340 条排除集，而非当前 1,350 条版本。108 个基因的另一次试算不属于本批结果。"},
 {id:"stage2",label:"蛋白互作证据整合",title:"蛋白互作证据整合",kicker:"实验互作记录 · 序列模型预测",summary:"A、B1 与 B2 使用 HI-union 与官方 RF2 筛选集合的并集；B3 使用异构体特异性 PLM/PPLM 预测。共保留 5,957 条来源记录，对应 5,954 个唯一蛋白配对。",version:"RF2-final80/90 是官方筛选集合，不能等同于任意配对原始分数超过 0.8/0.9。基因映射支持不能直接解释为某一异构体的实验验证。"},
 {id:"stage25",label:"同源与功能证据扩展",title:"同源与功能证据扩展",kicker:"候选蛋白与支持种子的证据关联",summary:"基于当前 A1/B1/B2 的 final80 未通过池进行同源功能扩展，独立列于附表，不计入直接候选–锚点配对。",version:"扩展采用当前 A1_WIDE、修正 B1 与重建 B2。A1_WIDE 不同于主候选路径 A 的 A3_FOCUSED；支持种子的互作结论不能转移至候选。"},
];
export const fmt=(n:number)=>n.toLocaleString("en-US");
export const pairKey=(p:Pair)=>[p.candidate,p.anchor].sort().join("|");
export function exportTSV(name:string,rows:Record<string,unknown>[]) {
 if(!rows.length)return;const cols=Object.keys(rows[0]);const clean=(v:unknown)=>String(v??"").replace(/[\t\r\n]/g," ");
 const text="\uFEFF"+[cols.join("\t"),...rows.map(r=>cols.map(c=>clean(r[c])).join("\t"))].join("\n");
 const url=URL.createObjectURL(new Blob([text],{type:"text/tab-separated-values;charset=utf-8"}));const a=document.createElement("a");a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
}
