"use client";

import {useEffect, useState, type ReactNode} from "react";
import {Atom, BookOpen, CircleDot, Dna, FileText, Fingerprint, GitMerge, Layers3, Link2, Microscope, Network, Orbit, Pause, Play, Search, ShieldCheck, Split, Waypoints} from "lucide-react";
import {Button} from "@/components/ui/button";
import {fmt, routeLabels, type Release, type Route} from "@/lib/release";
import {useBackgroundMotion} from '@/lib/use-background-motion';

type Props = {release:Release; search:ReactNode; query:string; go:(page:string,query?:string,route?:string)=>void; children:ReactNode};

// Intermediate counts are transcribed from the frozen method body, not new predictions.
// Final results and deduplicated counts always come from the protected release data.
const pathways = [
  {route:"A" as Route, id:"a", Icon:CircleDot, input:"3,418 条核定位注释蛋白", filter:"调控证据与功能新颖性分层", retained:2457, unit:"条聚焦候选", evidence:"HI-union ∪ RF2-final90", evidenceNote:"实验记录与官方预测集合并集"},
  {route:"B1" as Route, id:"b1", Icon:ScanProtein, input:"14,625 条非核注释代表序列", filter:"DeepLoc × NLSExplorer", retained:1350, unit:"条双模型支持候选", evidence:"HI-union ∪ RF2-final80", evidenceNote:"实验记录与官方预测集合并集"},
  {route:"B2" as Route, id:"b2", Icon:Waypoints, input:"13,275 条排除 B1 后的蛋白", filter:"四项网络指标 · 前 2%", retained:266, unit:"条网络邻近候选", evidence:"HI-union ∪ RF2-final80", evidenceNote:"实验记录与官方预测集合并集"},
  {route:"B3" as Route, id:"b3", Icon:Layers3, input:"保守分支 · 核定位与代表序列排除", filter:"异构体特异性核定位筛选", retained:497, unit:"条异构体候选", evidence:"PLM-interact > 0.99", evidenceSecond:"且 PPLM-PPI > 0.9", evidenceNote:"两模型联合判定"},
];

function ScanProtein({size=20}:{size?:number}) { return <Fingerprint size={size}/>; }

export default function ResearchOverview({release,search,query,go,children}:Props) {
  const [motion,setMotion]=useBackgroundMotion();
  const isoformCount=new Set(release.candidates.filter(c=>c.route==="B3").flatMap(c=>c.accessions)).size;
  const stats=[
    {label:"筛选路径记录",value:release.counts.candidate_rows,unit:"条",note:"保留同一基因的不同路径来源",Icon:Split,tone:"blue",page:"candidates"},
    {label:"候选基因",value:release.counts.unique_genes,unit:"个",note:"跨路径合并后的唯一基因",Icon:Dna,tone:"teal",page:"candidates"},
    {label:"候选蛋白条目",value:release.counts.unique_accessions,unit:"个",note:`含 ${isoformCount} 条可变异构体`,Icon:Atom,tone:"violet",page:"candidates"},
    {label:"唯一互作配对",value:release.counts.unique_pairs,unit:"对",note:`对应 ${fmt(release.counts.pair_rows)} 条来源记录`,Icon:Network,tone:"amber",page:"pairs"},
  ];
  return <>
    <div className="page-heading atlas-heading">
      <div><div className="eyebrow"><Orbit size={15}/>候选蛋白发现与证据分析</div><h1>表观遗传调控候选图谱</h1><p>整合核定位、调控网络与蛋白互作证据，追溯每一条候选的筛选依据。</p></div>
      <Button variant="outline" onClick={()=>go("downloads")}><BookOpen size={17}/>研究方法</Button>
    </div>
    <div className="overview-search">{search}<Button onClick={()=>go("candidates",query)}><Search size={17}/>检索候选</Button></div>
    <div className="stats-grid atlas-stats">
      {stats.map(({label,value,unit,note,Icon,tone,page})=><button className={`stat atlas-stat stat-${tone}`} key={label} onClick={()=>go(page)} aria-label={`查看${label}：${fmt(value)}${unit}`}><span className="stat-label">{label}<Icon size={25} strokeWidth={1.7}/></span><strong>{fmt(value)}<em>{unit}</em></strong><small>{note}</small></button>)}
    </div>

    <section className={`panel screening-panel ${motion?"motion-on":"motion-off"}`} aria-labelledby="screening-heading">
      <div className="chromatin-atmosphere" aria-hidden="true"><img src="/Protein_digger/images/chromatin-editorial.png" alt=""/></div>
      <div className="screening-header"><div><span className="eyebrow"><Waypoints size={15}/>核定位与调控关联的互补筛选</span><h2 id="screening-heading">多路径候选筛选</h2></div><div className="screening-controls"><span className="version-pill">四条互补筛选路径</span><button className="motion-toggle" aria-pressed={motion} onClick={()=>setMotion(v=>!v)} aria-label={motion?"暂停分子背景动效":"启用分子背景动效"} title={motion?"暂停背景动效":"启用背景动效"}>{motion?<Pause size={16}/>:<Play size={16}/>}<span>{motion?"暂停动效":"启用动效"}</span></button></div></div>

      <div className="screening-diagram">
        <div className="sequence-sources">
          <button className="sequence-source canonical-source" onClick={()=>go("stage0")}><Dna size={24}/><div><strong>人类审校代表序列</strong><span>核定位注释分区 · 排除已知调控锚点</span></div><b>20,416<small>条</small></b></button>
          <button className="sequence-source isoform-source" onClick={()=>go("b3")}><Layers3 size={22}/><div><strong>可变异构体序列</strong><span>独立序列集合</span></div><b>22,131<small>条</small></b></button>
        </div>
        <svg className="source-connectors" viewBox="0 0 1000 34" preserveAspectRatio="none" aria-hidden="true"><path d="M375 0 V12 H125 V34 M375 12 V34 M375 12 H625 V34"/><path className="isoform-connector" d="M875 0 V34"/><circle cx="375" cy="12" r="3"/></svg>
        <div className="pathway-grid">
          {pathways.map(p=><article key={p.route} className={`pathway pathway-${p.route}`}>
            <button className="pathway-method" onClick={()=>go(p.id)} aria-label={`查看${routeLabels[p.route]}方法`}>
              <div className="pathway-topline"><span className={`route-tag route-${p.route}`}>路径 {p.route}</span><p.Icon size={22}/></div>
              <h3>{routeLabels[p.route]}</h3><p className="pathway-input">{p.input}</p>
              <div className="pathway-filter">{p.filter}</div>
              <div className="pathway-intermediate"><b>{fmt(p.retained)}</b><span>{p.unit}</span></div>
              <span className="method-affordance"><BookOpen size={13}/>查看筛选方法</span>
            </button>
            <div className="evidence-connector" aria-hidden="true"/>
            <div className="pathway-evidence"><span className="evidence-caption"><Link2 size={14}/>互作证据筛选</span><strong>{p.evidence}{p.evidenceSecond&&<><br/>{p.evidenceSecond}</>}</strong><small>{p.evidenceNote}</small></div>
            <button className="pathway-result" onClick={()=>go("candidates","",p.route)} aria-label={`查看路径${p.route}的${release.counts.route_rows[p.route]}条筛选结果`}><span>筛选结果</span><b>{fmt(release.counts.route_rows[p.route])}<small>{p.route==="B3"?"个基因":"条记录"}</small></b><span className="result-sub">{p.route==="B3"?`${isoformCount} 条可变异构体`:`${fmt(release.counts.route_pair_rows[p.route])} 条配对来源记录`}</span></button>
          </article>)}
        </div>
        <svg className="merge-connectors" viewBox="0 0 1000 32" preserveAspectRatio="none" aria-hidden="true"><path d="M125 0 V15 H500 V32 M375 0 V15 M625 0 V15 M875 0 V15 H500"/><circle cx="500" cy="15" r="3"/></svg>
        <button className="candidate-convergence" onClick={()=>go("candidates")}><span className="convergence-label"><GitMerge size={23}/><span><strong>候选集合整合</strong><small>跨路径去重，保留全部证据来源</small></span></span><span className="convergence-count"><b>{fmt(release.counts.unique_genes)}</b> 个基因<span className="convergence-divider"/> {fmt(release.counts.candidate_rows)} 条路径记录</span></button>
      </div>
      <div className="screening-footnote"><ShieldCheck size={14}/><span>B3 采用保守分析分支，输入版本详见方法说明。分子背景为示意，不代表预测结构。</span></div>
    </section>

    <section className="evidence-reading" aria-label="互作与扩展证据的解释">
      <button onClick={()=>go("stage2")}><span className="evidence-icon experiment-icon"><Microscope size={24}/></span><div><h3>实验互作记录</h3><p>HI-union 的实验记录，经基因与蛋白标识映射关联候选。</p></div></button>
      <button onClick={()=>go("stage2")}><span className="evidence-icon prediction-icon"><Network size={24}/></span><div><h3>模型预测支持</h3><p>RF2、PLM 与 PPLM 提供互作预测，不等同于实验直接结合。</p></div></button>
      <button onClick={()=>go("stage25")}><span className="evidence-icon extension-icon"><GitMerge size={24}/></span><div><h3>同源与功能证据扩展 <b>{release.counts.appendix_rows}</b></h3><p>独立候选附表，不计入主集合，也不继承种子的直接互作结论。</p></div></button>
    </section>
    <section className="panel overview-candidates"><div className="panel-heading"><h2><FileText size={19}/>候选蛋白浏览</h2><Button variant="ghost" onClick={()=>go("candidates",query)}>查看完整候选集</Button></div>{children}</section>
  </>;
}
