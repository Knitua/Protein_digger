"use client";
import {Dna,Split,Network,GitMerge,MapPin,Pause,Play,CircleDot,Fingerprint,Waypoints,Layers3,Info,type LucideIcon} from 'lucide-react';
import {workflowBranches} from '@/lib/workflow';
import {useBackgroundMotion} from '@/lib/use-background-motion';

function Connector({kind='straight',active=false}:{kind?:'straight'|'split'|'merge';active?:boolean}){
 const paths=kind==='split'
  ? ['M 1 94 H 20 Q 32 94 32 82 V 56 Q 32 44 44 44 H 60','M 32 94 V 132 Q 32 144 44 144 H 60']
  : kind==='merge'?['M 1 44 H 20 Q 32 44 32 56 V 94','M 1 144 H 20 Q 32 144 32 132 V 94','M 32 94 H 60']:['M 1 94 H 60'];
 const ends=kind==='split'?[44,144]:[94];
 return <div className={`ribbon-connector ribbon-connector-${kind} ${active?'is-active':''}`} aria-hidden="true">
  <svg className="ribbon-connector-wide" viewBox="0 0 64 188" preserveAspectRatio="none">
   {paths.map((d,i)=><path key={i} d={d}/>)}
   {ends.map(y=><path key={y} className="ribbon-arrowhead" d={`M 55 ${y-4} L 61 ${y} L 55 ${y+4}`}/>)}
  </svg>
  {kind==='merge'&&<span className="ribbon-mobile-source">互作候选与同源扩展共同收录</span>}
  <svg className="ribbon-connector-narrow" viewBox="0 0 40 32"><path d="M20 1 V28 M15 23 L20 28 L25 23"/></svg>
 </div>;
}

export default function WorkflowMap({page,go}:{page:string;go:(p:string)=>void}){
 const [motion,setMotion]=useBackgroundMotion();
 const icons=[CircleDot,Fingerprint,Waypoints,Layers3];
 const inBranch=workflowBranches.some(b=>b.id===page);
 const node=(id:string,title:string,stage:string,detail:string,Icon:LucideIcon)=><button
  data-flow-node={id} className={`ribbon-node ${id==='stage25'?'ribbon-extension':''} ${page===id?'is-current':''}`}
  aria-current={page===id?'step':undefined} onClick={()=>go(id)}>
  <span className="ribbon-icon"><Icon size={21} strokeWidth={1.65}/></span>
  <span className="ribbon-node-copy"><small>{stage}</small><strong>{title}</strong>
   <span className={`ribbon-node-note ${page===id?'current-note':''}`}>{page===id?<><i/>当前查看</>:detail}</span>
  </span>
 </button>;
 return <section className={`panel workflow-ribbon ${motion?'motion-on':'motion-off'}`} aria-label="研究流程定位">
  <div className="ribbon-atmosphere" aria-hidden="true">
   <div className="chromatin-atmosphere"><img src="/Protein_digger/images/chromatin-editorial.png" alt=""/></div>
   <div className="ribbon-light ribbon-light-blue"/><div className="ribbon-light ribbon-light-teal"/>
   <svg className="ribbon-molecular-field" viewBox="0 0 1200 300" preserveAspectRatio="xMidYMid slice">
    <path d="M-80 280 C60 90 150 335 310 150 S540 265 705 85 S1015 210 1280 -15"/>
    <path d="M-70 250 C80 365 170 110 325 220 S520 15 725 150 S1020 -60 1260 90"/>
    {[110,275,440,605,770,935,1100].map((x,i)=><g key={x} transform={`translate(${x} ${175+Math.sin(i*1.8)*75}) rotate(-25)`}><ellipse rx="28" ry="14"/><ellipse rx="19" ry="24"/><circle r="8"/></g>)}
   </svg>
  </div>
  <header className="ribbon-heading"><span className="ribbon-heading-icon"><MapPin size={17}/></span><h2>研究流程</h2><span className="ribbon-heading-hint">从蛋白组到候选集合 · 点击步骤切换</span><button className="ribbon-motion" aria-pressed={motion} aria-label={motion?'暂停分子背景动效':'启用分子背景动效'} onClick={()=>setMotion(v=>!v)}>{motion?<Pause size={14}/>:<Play size={14}/>}</button></header>
  <div className="ribbon-grid">
   {node('reference','蛋白组分区','Stage0','代表序列与调控锚点',Dna)}
   <Connector active={page==='reference'||inBranch}/>
   <div className="ribbon-paths">{workflowBranches.map((b,i)=>{const Icon=icons[i];return <button key={b.id} data-flow-node={b.id} className={`ribbon-path ribbon-path-${b.route} ${page===b.id?'is-current':''}`} aria-current={page===b.id?'step':undefined} onClick={()=>go(b.id)}><span className="ribbon-path-icon"><Icon size={18} strokeWidth={1.65}/></span><span><small>路径 {b.route}</small><strong>{b.label}</strong>{page===b.id&&<em>当前查看</em>}</span></button>})}<span className="ribbon-path-caption">Stage1 · 四条互补路径</span></div>
   <Connector kind="split" active={inBranch||page==='pairs'||page==='stage25'}/>
   <div className="ribbon-evidence">
    {node('pairs','互作证据筛选','Stage2','支持直接配对的候选',Network)}
    {node('stage25','同源功能扩展','Stage2','独立扩展附表',Split)}
   </div>
   <Connector kind="merge" active={page==='pairs'||page==='stage25'||page==='candidates'}/>
   {node('candidates','候选集合整合','跨路径归并','主集合与独立扩展附表',GitMerge)}
  </div>
  <footer className="ribbon-boundary"><Info size={13}/><span>互作候选与同源扩展分别收录于主集合和独立附表；同源支持不等同于直接互作证据。</span></footer>
 </section>;
}
