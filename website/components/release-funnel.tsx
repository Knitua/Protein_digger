"use client";
import {Filter} from 'lucide-react';
import {Release,Route,fmt,routeLabels} from '@/lib/release';
import {useAtlas,RecordRow} from '@/lib/atlas';
import ReaderNote from '@/components/reader-note';

export default function ReleaseFunnel({release,route,onRoute}:{release:Release;route:string;onRoute:(r:string)=>void}){
 const {data}=useAtlas<RecordRow>('workflow_summary');
 return <section className="panel local-funnel release-summary">
  <div className="local-funnel-heading"><h2>各路径互作筛选</h2></div>
  <div className="release-route-grid">{(['A','B1','B2','B3'] as Route[]).map(r=>
   <button key={r} className={`release-route-card route-${r} ${route===r?'selected':''}`} aria-pressed={route===r} aria-label={`筛选路径 ${r}：${routeLabels[r]}`} onClick={()=>onRoute(route===r?'all':r)}>
    <div className="release-route-heading"><span className="release-route-badge">{r}</span><h3>{routeLabels[r]}</h3><Filter size={14}/></div>
    <div className="release-route-counts"><div><span>入选{r==='B3'?'基因':'蛋白'}</span><strong>{fmt(release.counts.route_rows[r])}<small>{r==='B3'?'个':'条'}</small></strong></div><div className="release-route-input"><span>输入{r==='B3'?'异构体':'蛋白'}</span><b>{data?fmt(data.stage1[r]):'…'}</b></div></div>
    <div className="release-route-rule">{r==='B3'?<><span>PLM-interact &gt; 0.99</span><i>且</i><span>PPLM-PPI &gt; 0.9</span></>:<><span>HI-union</span><i>∪</i><span>{r==='A'?'RF2-final90':'RF2-final80'}</span></>}</div>
    <div className="release-route-footer">{r==='B3'&&<span><b>{fmt(new Set(release.candidates.filter(c=>c.route==='B3').flatMap(c=>c.accessions)).size)}</b> 条异构体</span>}<span><b>{fmt(release.counts.route_pair_rows[r])}</b> 条配对来源</span></div>
   </button>)}</div>
  <ReaderNote title="筛选口径与证据边界"><p>各路径保留原始证据来源，跨路径共有 {fmt(release.counts.unique_pairs)} 个唯一蛋白配对。HI-union 与 RF2-final80/90 采用并集；B3 的两项模型条件必须同时满足。不同模型分数不合成为统一互作强度。B3 分别统计基因与异构体，不能直接相加。</p></ReaderNote>
 </section>;
}
