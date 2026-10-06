"use client";
import '@/app/nucleus-model.css';
import {useEffect,useState} from 'react';
import {ArrowRight,Download,ExternalLink,Focus,Layers3} from 'lucide-react';
import {Button} from '@/components/ui/button';
import NucleusArchitecture from '@/components/nucleus-architecture';
import {benchmarkBase,validateNucleusBenchmark,type NucleusBenchmark} from '@/lib/nucleus-benchmark';

const modules=[
 {id:'encoder',name:'序列表征',text:'ESM2-3B 编码器冻结，输出第 36 层的逐残基特征；长序列保留两端各 511 个残基。'},
 {id:'general',name:'通用分支',text:'11 个可学习标签查询分别汇聚序列信息，输出十类亚细胞定位与一个辅助膜标签。'},
 {id:'specialist',name:'专用分支',text:'三个查询分别接受 NLS、NES 和非核信号约束，经共享专家映射组合为核定位专用输出。注意力不等同于已验证的定位位点。'},
 {id:'fusion',name:'门控融合',text:'可学习门控将专用输出叠加至通用核定位输出，其余标签不变；实际组合为 NLS 专家 + NES 专家 − 非核专家。'},
];

export default function NucleusModelWorkspace({onB1}:{onB1:()=>void}){
 const [data,setData]=useState<NucleusBenchmark|null>(null),[error,setError]=useState(''),[attempt,setAttempt]=useState(0);
 const [part,setPart]=useState('specialist'),[axis,setAxis]=useState<'detail'|'zero'>('detail');
 useEffect(()=>{
  const controller=new AbortController();setError('');
  fetch(`${benchmarkBase}/benchmark-v1.json`,{cache:'no-store',signal:controller.signal})
   .then(r=>{if(!r.ok)throw Error('模型数据暂时无法载入。');return r.json()})
   .then(validateNucleusBenchmark).then(setData)
   .catch(e=>{if(e.name!=='AbortError')setError(e.message)});
  return()=>controller.abort();
 },[attempt]);
 if(!data)return <section className="panel nm-loading" role="status">{error?<><p>{error}</p><Button variant="outline" onClick={()=>setAttempt(x=>x+1)}>重新载入</Button></>:'正在载入模型数据…'}</section>;
 const selectedPart=modules.find(m=>m.id===part)!;
 const chart=`${benchmarkBase}/figures/benchmark-bars-${axis}-v2.svg`;
 return <div className="nucleus-model-workspace">
  <header className="nm-identity">
   <div><span className="nm-kicker"><Focus size={16}/>核定位专用模型</span><h2>Nucleus Specialist</h2><p>共享序列表征，双分支学习，核定位门控融合。</p></div>
   <button onClick={onB1}>尚未用于当前 B1 筛选<ArrowRight size={15}/></button>
  </header>
  <section className="nm-panel" aria-labelledby="nm-architecture-title">
   <div className="nm-section-heading"><h2 id="nm-architecture-title"><Layers3 size={20}/>模型架构</h2><span>冻结编码器 · 可训练预测头</span></div>
   <div className="nm-architecture-scroll" tabIndex={0} role="region" aria-label="模型架构示意，窄屏可横向滚动">
    <NucleusArchitecture selected={part} onSelect={setPart}/>
   </div>
   <div className="nm-module-strip">
    <div className="nm-module-tabs" role="group" aria-label="查看模型模块">{modules.map(m=><button key={m.id} aria-pressed={part===m.id} onClick={()=>setPart(m.id)}>{m.name}</button>)}</div>
    <p aria-live="polite"><strong>{selectedPart.name}</strong>{selectedPart.text}</p>
   </div>
  </section>
  <section className="nm-panel nm-benchmark" aria-labelledby="nm-benchmark-title">
   <div className="nm-section-heading"><h2 id="nm-benchmark-title">与 DeepLoc 2.0 比较</h2><div className="nm-chart-actions"><div className="nm-axis-toggle" role="group" aria-label="柱状图纵轴范围"><button aria-pressed={axis==='detail'} onClick={()=>setAxis('detail')}>局部纵轴</button><button aria-pressed={axis==='zero'} onClick={()=>setAxis('zero')}>完整纵轴</button></div><a href={chart} download aria-label="下载当前柱状图 SVG"><Download size={17}/></a></div></div>
   <div className="nm-chart-scroll" tabIndex={0} role="region" aria-label="性能比较柱状图，窄屏可横向滚动"><img src={chart} width="1440" height="510" alt={`DeepLoc 2.0 与 Nucleus Specialist 的七项指标，均值、标准差与五折散点；纵轴 ${axis==='detail'?'0.60–1.00，非零起点':'0–1，零基线'}。`}/></div>
   <div className="nm-chart-note"><span>选定指标展示；完整结果见下载资料。模型差异不等同于统计显著性。</span><a href={`${benchmarkBase}/metrics.tsv`} download>完整指标数据<ArrowRight size={14}/></a></div>
  </section>
  <footer className="nm-resources"><span>方法与数据</span><a href={`${benchmarkBase}/method.md`} download><Download size={15}/>模型方法</a><a href={`${benchmarkBase}/benchmark-v1.json`} download><Download size={15}/>原始指标与来源</a><a href="https://github.com/Knitua/Protein_digger/tree/d91687d/benchmarks/deeploc" target="_blank" rel="noreferrer"><ExternalLink size={15}/>代码与基准</a></footer>
 </div>;
}
