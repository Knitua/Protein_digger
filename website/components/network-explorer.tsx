"use client";
import {useEffect,useMemo,useRef,useState} from 'react';
import {Search,LocateFixed,RotateCcw,ChevronLeft,Plus} from 'lucide-react';
import {useAtlas,RecordRow} from '@/lib/atlas';
import {fmt} from '@/lib/release';
import {Button} from '@/components/ui/button';
import {Input} from '@/components/ui/input';
import {Tabs,TabsList,TabsTrigger} from '@/components/ui/tabs';
import AtlasGraph,{AtlasNode,AtlasEdge,endpoint} from '@/components/atlas-graph';
import ReaderNote from '@/components/reader-note';
import NetworkOverview from '@/components/network-overview';
type Result={nodes:AtlasNode[];links:AtlasEdge[];availableNodes:number;availableEdges:number;covered:boolean;neighborCount:number;pathStatus:string;path:string[];limited:boolean};
export default function NetworkExplorer({focus,onFocus,onProtein,query,onQuery}:{focus:string;onFocus:(id:string)=>void;onProtein:(id:string)=>void;query:string;onQuery:(s:string)=>void}){
 const {data:reference}=useAtlas<RecordRow[]>('reference');const [mode,setMode]=useState('annotation'),[scope,setScope]=useState('one'),[limit,setLimit]=useState(1),[result,setResult]=useState<Result|null>(null),[error,setError]=useState(''),[loading,setLoading]=useState(false),[edge,setEdge]=useState<AtlasEdge|null>(null),[trail,setTrail]=useState<string[]>([]),[target,setTarget]=useState(''),[term,setTerm]=useState('');
 const studio=useRef<HTMLElement>(null);
 const worker=useRef<Worker|null>(null),request=useRef(0);const meta=useMemo(()=>new Map(reference?.map(r=>[r.accession,r])||[]),[reference]);
 useEffect(()=>{const w=new Worker('/Protein_digger/network-worker.js',{type:'module'});worker.current=w;w.onmessage=({data})=>{if(data.requestId!==request.current)return;setLoading(false);setError(data.error||'');setResult(data.result||null)};w.onerror=()=>{setLoading(false);setError('网络索引载入失败，请重新载入页面；数据库仍可使用。')};return()=>{worker.current=null;w.terminate()}},[]);
 useEffect(()=>{setLimit(1);setEdge(null);setTerm('');setTarget('')},[focus,mode,scope]);
 useEffect(()=>{if(!focus){request.current++;setResult(null);setLoading(false);return;}setLoading(true);setResult(null);worker.current?.postMessage({requestId:++request.current,query:{mode,focus,scope,nodeLimit:150*limit,edgeLimit:400*limit,target}})},[focus,mode,scope,limit,target]);
 const choose=(id:string)=>{if(id===focus)return;setTrail(t=>[...t,focus].filter(Boolean));onFocus(id)};
 const options=useMemo(()=>{const q=query.trim().toLowerCase();if(!q)return [];return (reference||[]).filter(r=>`${r.gene} ${r.accession} ${r.name}`.toLowerCase().includes(q)).slice(0,10)},[reference,query]);
 // Align canvases after their responsive toolbars wrap; not a data transform.
 useEffect(()=>{
  const el=studio.current;if(!el)return;
  let frame=0;
  const align=()=>{cancelAnimationFrame(frame);frame=requestAnimationFrame(()=>{
   const panes=[el.querySelector<HTMLElement>('.density-canvas'),el.querySelector<HTMLElement>('.network-center-pane .graph-canvas')];
   if(panes.some(p=>!p))return;
   const boxes=panes as HTMLElement[],wide=matchMedia('(min-width:1121px)').matches;
   const tops=boxes.map(p=>p.getBoundingClientRect().top-(parseFloat(getComputedStyle(p).marginTop)||0)),max=Math.max(...tops);
   boxes.forEach((p,i)=>{const next=wide?Math.max(0,max-tops[i]):0;if(Math.abs((parseFloat(p.style.marginTop)||0)-next)>.5)p.style.marginTop=next+'px'});
  })};
  const observer=new ResizeObserver(align);
  el.querySelectorAll('.network-split,.atlas-graph-header,.density-toolbar,.network-controls,.network-coverage').forEach(n=>observer.observe(n));
  align();return()=>{cancelAnimationFrame(frame);observer.disconnect()};
 },[result,mode,scope]);
 const targetRecord=target?meta.get(target):null;
 return <section ref={studio} className="panel network-studio regulatory-studio"><div className="network-studio-controls"><div className="network-source-switch"><Tabs value={mode} onValueChange={setMode}><TabsList><TabsTrigger value="annotation">功能注释关联网络</TabsTrigger><TabsTrigger value="omni">OmniPath 有向调控网络</TabsTrigger></TabsList></Tabs></div>
 <section className="network-focus-bar" id="network-center"><div><LocateFixed size={23}/><span><small>中心蛋白</small><strong>{meta.get(focus)?.gene||focus||'请选择蛋白'}</strong><code>{focus}</code></span></div><div className="network-focus-actions"><Button variant="outline" disabled={!trail.length} onClick={()=>{const id=trail.at(-1)!;setTrail(t=>t.slice(0,-1));onFocus(id)}}><ChevronLeft size={15}/>上一蛋白</Button><Button variant="outline" disabled={!focus} onClick={()=>onProtein(focus)}>查看蛋白详情</Button></div></section>
 <div className="network-search"><div className="search-input"><Search size={17}/><Input aria-label="检索网络中心蛋白" placeholder="输入基因或UniProt编号，选择中心蛋白" value={query} onChange={e=>onQuery(e.target.value)}/></div>{query&&<div className="network-search-results">{options.length?options.map(r=><button key={r.accession} className={r.accession===focus?'selected':''} onClick={()=>choose(r.accession)}><strong>{r.gene||'无基因符号'}</strong><span>{r.accession}</span></button>):<span>没有匹配的参考蛋白；不会保留上一次蛋白的子图。</span>}</div>}</div>
 </div><div className="network-split"><div className="network-overview-pane"><NetworkOverview mode={mode} focus={focus} onSelect={choose} labels={meta}/></div><div className="network-center-pane"><div className="network-controls"><Tabs value={scope} onValueChange={setScope}><TabsList><TabsTrigger value="one">一跳邻域</TabsTrigger><TabsTrigger value="two">两跳邻域</TabsTrigger><TabsTrigger value="path">通向锚点的最短路径</TabsTrigger></TabsList></Tabs>{scope==='path'&&<label className="path-target">锚点编号（可选）<Input aria-label="目标锚点编号" placeholder="留空：最近可达锚点" value={target} onChange={e=>setTarget(e.target.value.trim().toUpperCase())}/></label>}</div>
 {scope==='path'&&target&&!targetRecord?.anchor?<p className="notice">请输入工作锚点的有效UniProt编号；不会把普通蛋白作为锚点。</p>:null}
 {loading?<div className="notice">正在查询真实网络关系，下方数据库仍可操作。</div>:null}{error&&<p className="notice">{error}</p>}
 {result&&!(scope==='path'&&target&&!targetRecord?.anchor)&&<><div className="network-coverage"><span>已展示 {fmt(result.nodes.length)} / {fmt(result.availableNodes)} 个节点 · {fmt(result.links.length)} / {fmt(result.availableEdges)} 条关系</span>{result.limited&&<Button variant="outline" onClick={()=>setLimit(n=>n+1)}><Plus size={15}/>继续展开</Button>}{limit>1&&<Button variant="ghost" onClick={()=>setLimit(1)}><RotateCcw size={14}/>恢复初始范围</Button>}</div>
 {!result.covered?<p className="notice">当前资源未覆盖这个蛋白，保留中心节点供辨认；这不等同于没有生物学联系。</p>:!result.neighborCount?<p className="notice">当前资源有该节点，但未记录相邻关系。</p>:null}
 {result.pathStatus==='NO_REACHABLE_ANCHOR'&&<p className="notice">按当前网络及其方向未找到可达锚点，不解释为不存在调控关联。</p>}
 <AtlasGraph embedded showCounts={false} height={430} nodes={result.nodes} links={result.links} selected={focus} centerSelected title={mode==='annotation'?'蛋白中心 · 功能注释关系':'蛋白中心 · 有向调控关系'} directed={mode==='omni'} allow3D={false} onNode={id=>{if(id.startsWith('term:')){setTerm(id);setEdge(null)}else {choose(id);onProtein(id)}}} onEdge={e=>{setEdge(e);setTerm('')}}/>
 {scope==='path'&&result.pathStatus==='FOUND'&&<div className="network-path"><strong>真实边组成的最短路径</strong><p>{result.path.map(id=>meta.get(id)?.gene||result.nodes.find(n=>n.id===id)?.label||id).join(mode==='omni'?' → ':' — ')}</p></div>}
 
 </>}
 {!focus&&!loading&&<p className="notice">当前检索范围没有中心蛋白，请重新检索。</p>}
 </div></div><ReaderNote title="网络含义与显示说明"><p>注释网络的一跳邻域表示蛋白与功能术语的关联，两跳展开共享术语的蛋白与锚点；共同注释不等于直接互作。OmniPath 连线保留资源中的调控方向。</p><p>热图仅表示网络布局网格内的节点数量，采用对数色阶，每级使用全图固定色域，不随平移变化。锚点是蛋白子集，不重复相加。颜色与布局不表示功能活性、富集程度、互作强度或生物学空间位置。</p><p>中心图的展示上限不改变完整邻域、筛选排名或已选择的最短路径；需要时可继续展开。</p></ReaderNote><div className="network-evidence-details">
 {term&&<div className="notice"><strong>功能术语</strong><p>{result?.nodes.find(n=>n.id===term)?.label||term}</p><p>节点表示功能注释，不自动代表与中心蛋白直接相连；点击真实连线查看关联来源。</p></div>}{edge&&<div className="notice"><div><strong>{edge.kind==='annotation'?'功能注释关联':'有向调控记录'}</strong><p>{edge.label||'来源未提供'}</p><p>{meta.get(endpoint(edge.source))?.gene||endpoint(edge.source)} {edge.kind==='omnipath'?'→':'—'} {meta.get(endpoint(edge.target))?.gene||result?.nodes.find(n=>n.id===endpoint(edge.target))?.label||endpoint(edge.target)}</p></div></div>}
 </div></section>;
}
