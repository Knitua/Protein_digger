"use client";
import {ArrowRight,Filter,Layers3,Database,GitMerge} from 'lucide-react';
import {RecordRow,regulatoryLabels,tierLabels,useAtlas} from '@/lib/atlas';
import ReaderNote from '@/components/reader-note';
import {fmt} from '@/lib/release';
type Props={page:string;data:RecordRow[];group:string;secondary:string;onSelect:(g:string,s?:string)=>void};
type AnchorComposition={sources:{id:string;name:string;count:number;description:string}[];lambertAnimalUnion:number;lambertAnimalOverlap:number;epifactorsOverlap:number;union:number;animalStrict:number;animalExtension:number};
export default function StageFunnel({page,data,group,secondary,onSelect}:Props){
 const {data:anchors,error:anchorError}=useAtlas<AnchorComposition>(page==='reference'?'anchor_composition':null);
 const n=(f:(r:RecordRow)=>boolean)=>data.filter(f).length;
 const item=(id:string,label:string,value:number,note='',sub='all')=>{
  const selected=group===id&&(sub==='all'||secondary===sub);
  return <button className={`funnel-node ${selected?'selected':''}`} aria-pressed={selected} onClick={()=>onSelect(id,sub)}><span>{label}</span><strong>{fmt(value)}<small>条</small></strong>{note&&<p>{note}</p>}<span className="funnel-action"><Filter size={12}/><span>查看记录</span></span></button>;
 };
 const arrow=<div className="compact-flow-arrow" aria-hidden="true"><ArrowRight size={19}/></div>;
 const focused=n(r=>r.regulatory_evidence_status==='KNOWN_COMPLEX_ASSOCIATED'||r.regulatory_evidence_status==='NO_KNOWN_REGULATORY_EVIDENCE'&&['UNDER_CHARACTERIZED','MIXED_MOONLIGHTING'].includes(r.final_tier));
 return <section className={`panel local-funnel compact-funnel funnel-${page}`} aria-label="当前步骤筛选流程">
  <div className="local-funnel-heading"><Layers3 size={18}/><h2>{({reference:'核定位注释分区与锚点排除',a:'路径 A · 调控证据与功能分层',b1:'路径 B1 · 双模型核定位筛选',b2:'路径 B2 · 网络邻近性筛选',b3:'路径 B3 · 异构体独立评价'} as Record<string,string>)[page]}</h2></div>
  {page==='reference'&&<>
   <div className="reference-sources">
    <div className="proteome-source"><h3><Database size={15}/>参考蛋白组</h3>{item('all','Swiss-Prot 人类审校代表序列',data.length,'核定位注释须适用于代表序列')}</div>
    <div className="anchor-composition"><h3><GitMerge size={15}/>调控锚点数据集 <small>筛选后蛋白数</small></h3>
     {anchors?<><div className="anchor-source-grid">{anchors.sources.map(s=><div key={s.id} className={`anchor-source anchor-${s.id}`}><span>{s.name}</span><strong>{fmt(s.count)}<small>条</small></strong><p>{s.description}</p></div>)}<button className={`anchor-union ${group==='anchors'?'selected':''}`} aria-pressed={group==='anchors'} onClick={()=>onSelect('anchors')}><span>三来源去重并集</span><strong>{fmt(n(r=>r.anchor))}<small>条</small></strong><p><Filter size={11}/>查看全部锚点</p></button></div><ReaderNote title="来源组成与去重"><p>Lambert 与 AnimalTFDB 接受集重叠 {fmt(anchors.lambertAnimalOverlap)} 条，合并为 {fmt(anchors.lambertAnimalUnion)} 条；再与 EpiFactors 合并时重叠 {fmt(anchors.epifactorsOverlap)} 条。AnimalTFDB 接受集含 {anchors.animalStrict} 条严格 TF 支持、{anchors.animalExtension} 条调控扩展，非全库数量。</p></ReaderNote></>:<p className="anchor-overlap">{anchorError||'正在载入来源组成…'}</p>}
    </div>
   </div>
   <div className="partition-lanes">{(['A','B'] as const).map(p=><div className={`partition-lane lane-${p}`} key={p}><span className="lane-label">{p}</span>{item(p,p==='A'?'具有适用核定位注释':'缺少适用核定位注释',n(r=>r.partition===p))}{arrow}{item('excluded-'+p,'按锚点规则排除',n(r=>r.partition===p&&r.excluded))}{arrow}{item(p==='A'?'a-input':'b-input',p==='A'?'A 路线输入':'B1 / B2 共同输入',n(r=>r.partition===p&&!r.excluded))}</div>)}</div>
   <ReaderNote title="分区与排除规则"><p>A、B 全集互斥且覆盖全部代表序列。锚点按蛋白编号及独立基因符号排除，因此排除记录数不等于锚点并集数。可变异构体独立评价。</p></ReaderNote>
  </>}
  {page==='a'&&<>
   <div className="funnel-input-strip">{item('all','已注释核定位、非锚点蛋白',data.length)}<p>先按调控证据分为四类，再对无已知调控证据的蛋白作功能分层。</p></div>
   <div className="funnel-four compact-evidence-grid">{Object.entries(regulatoryLabels).map(([id,label],i)=><div key={id} className={`evidence-class evidence-class-${i}`}>{item(id,label,n(r=>r.regulatory_evidence_status===id))}</div>)}</div>
   <div className="compact-tier-section"><div><span>无已知调控证据</span><strong>二级功能分层</strong></div><div className="funnel-tier-list">{Object.entries(tierLabels).filter(([t])=>n(r=>r.regulatory_evidence_status==='NO_KNOWN_REGULATORY_EVIDENCE'&&r.final_tier===t)>0).map(([t,label])=><button key={t} aria-pressed={group==='NO_KNOWN_REGULATORY_EVIDENCE'&&secondary===t} className={group==='NO_KNOWN_REGULATORY_EVIDENCE'&&secondary===t?'selected':''} onClick={()=>onSelect('NO_KNOWN_REGULATORY_EVIDENCE',t)}><span>{label}</span><b>{fmt(n(r=>r.regulatory_evidence_status==='NO_KNOWN_REGULATORY_EVIDENCE'&&r.final_tier===t))}</b></button>)}</div></div>
   <div className="funnel-output-strip"><p>保留调控关联、研究不足及混合功能三类。</p>{arrow}{item('focused','聚焦候选集合',focused)}</div>
  </>}
  {page==='b1'&&<><div className="compact-screening-chain">{item('all','无适用核注释、非锚点代表序列',data.length)}{arrow}<div className="compact-model-branches">{item('deeploc','DeepLoc 核定位支持',n(r=>r.deeploc_gt_0p5),'核定位分数 > 0.5')}{item('nls','NLSExplorer 窗口支持',n(r=>r.nlsexplorer_gt_0p5),'最高分窗口分数 > 0.5')}</div>{arrow}{item('positive','双模型共同支持',n(r=>r.both_gt_0p5),'两项条件同时满足')}</div><ReaderNote title="定位证据说明"><p>两条评分支路取交集，不相加；模型最高分窗口不等同于经过实验验证的核定位信号。</p></ReaderNote></>}
  {page==='b2'&&<><div className="compact-screening-chain">{item('all','B 路线完整网络特征',data.length)}{arrow}<div className="compact-model-branches">{item('excluded','B1 已覆盖，排除排序',n(r=>r.b1Excluded))}{item('eligible','正式网络排序背景',n(r=>!r.b1Excluded))}</div>{arrow}{item('selected','调控网络邻近候选',n(r=>r.selected),'四项百分位等权聚合 · 前 2%')}</div><ReaderNote title="网络指标说明"><p>四项指标：注释图距离、共享种子 IDF、有向距离、可达种子数。采用冻结排名；网络关联不代表物理互作。</p></ReaderNote></>}
  {page==='b3'&&<><div className="compact-isoform-chain">{item('all','真实可变异构体全集',data.length)}{arrow}{item('parent','保守分支母本范围',n(r=>r.parentEligible),'母本代表序列缺少适用核注释')}{arrow}{item('selected','保守分支异构体候选',n(r=>r.selected),'定位支持及基因、锚点排除')}{arrow}{item('delivery','联合互作条件通过',n(r=>r.delivery.includes('B3')),'PLM-interact > 0.99 且 PPLM-PPI > 0.9')}</div><ReaderNote title="保守分支说明"><p>异构体定位证据／双模型支持后，排除历史 B1 已覆盖基因与调控锚点。保留该分支的历史 B1 排除集合，未评分不作阴性处理。</p></ReaderNote></>}
 </section>;
}
