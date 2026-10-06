"use client";
import {useMemo,useState} from 'react';
import {ColumnDef} from '@tanstack/react-table';
import {Search,Download,SlidersHorizontal} from 'lucide-react';
import {Input} from '@/components/ui/input';
import {Select,SelectContent,SelectItem,SelectTrigger,SelectValue} from '@/components/ui/select';
import {Button} from '@/components/ui/button';
import DataTable from '@/components/data-table';
import HomologyRelations from '@/components/homology-relations';
import HomologyEvidence from '@/components/homology-evidence';
import AtlasMethods from '@/components/atlas-methods';
import {useAtlas,RecordRow,display} from '@/lib/atlas';
import {homologyNumber,homologyPercent} from '@/lib/homology-display';
import {Doc,exportTSV,fmt} from '@/lib/release';

type ExtensionSummary={screened_canonical_universe:number;primary_final80_positive_seeds:number;primary_final80_failure_pool:number;input_routes:Record<string,{screened:number;final80_positive:number;final80_failure:number}>};

export default function ExtensionAtlas({query,onQuery,onProtein,docs}:{query:string;onQuery:(s:string)=>void;onProtein:(s:string)=>void;docs:Doc[]}){
 const {data,error}=useAtlas<{records:RecordRow[];edges:RecordRow[];summary:ExtensionSummary}>('extension');
 const [route,setRoute]=useState('all'),[functional,setFunctional]=useState('all');
 const [selected,setSelected]=useState<RecordRow|null>(null),[fullMetrics,setFullMetrics]=useState(false);
 const rows=useMemo(()=>data?.records.filter(r=>
  (route==='all'||r.candidate_stage1_route===route)&&
  (functional==='all'||functional==='mf'&&r.shared_go_mf_direct_ids||functional==='reactome'&&r.shared_reactome_direct_ids||functional==='bpmf'&&r.functional_pass_reason.includes('bpmf'))&&
  `${r.candidate_gene} ${r.candidate_accession} ${r.best_seed_gene} ${r.best_seed_accession}`.toLowerCase().includes(query.toLowerCase())
 )||[],[data,query,route,functional]);
 const active=rows.find(r=>r.candidate_accession===selected?.candidate_accession)||[...rows].sort((a,b)=>a.candidate_gene.localeCompare(b.candidate_gene)||a.candidate_accession.localeCompare(b.candidate_accession))[0]||null;
 const cols:ColumnDef<RecordRow>[]=[
  {accessorKey:'candidate_gene',header:'扩展候选',cell:({row})=><button className="gene-link" onClick={()=>onProtein(row.original.candidate_accession)}>{row.original.candidate_gene}<small className="cell-sub">{row.original.candidate_accession}</small></button>},
  {accessorKey:'best_seed_gene',header:'支持种子',cell:({row})=><button className="text-link" onClick={()=>onProtein(row.original.best_seed_accession)}>{row.original.best_seed_gene}</button>},
  {id:'sequence_identity',accessorFn:r=>homologyNumber(r.sequence_identity),header:'序列一致性',cell:({getValue})=>homologyPercent(getValue())},
  ...(fullMetrics?(['candidate_coverage','seed_coverage'] as const).map((key,i)=>({id:key,accessorFn:(r:RecordRow)=>homologyNumber(r[key]),header:['候选覆盖率','种子覆盖率'][i],cell:({getValue}:any)=>homologyPercent(getValue())})):[]),
  {id:'domain_boundary_min_normalized_iou',accessorFn:r=>homologyNumber(r.domain_boundary_min_normalized_iou),header:'结构域边界 IoU',cell:({getValue})=>display(getValue())},
  {id:'evidence',header:'支持证据',cell:({row})=><button className="text-link" onClick={()=>{setSelected(row.original);document.querySelector('.homology-comparison')?.scrollIntoView({block:'center',behavior:window.matchMedia('(prefers-reduced-motion: reduce)').matches?'auto':'smooth'})}}>查看对照</button>}
 ];
 if(!data)return <div className="panel loading">{error||'正在载入同源功能扩展…'}</div>;
 return <>
  <section className="panel local-funnel"><div className="local-funnel-heading"><h2>独立候选—种子支持流程</h2></div><div className="release-funnel">
   <button onClick={()=>{setRoute('all');setFunctional('all');onQuery('')}}><span>当前已筛查蛋白</span><strong>{fmt(data.summary.screened_canonical_universe)} 条</strong><small>{Object.entries(data.summary.input_routes).map(([r,c])=>`${r} ${fmt(c.screened)}`).join(' · ')}</small></button>
   <button onClick={()=>setFunctional('all')}><span>final80 未通过池</span><strong>{fmt(data.summary.primary_final80_failure_pool)} 条</strong><small>与 {fmt(data.summary.primary_final80_positive_seeds)} 条阳性种子比较</small></button>
   <button onClick={()=>{setRoute('all');setFunctional('all');onQuery('')}}><span>独立扩展候选</span><strong>{data.records.length} 条</strong><small>序列、结构域及严格功能证据联合支持</small></button>
  </div><p className="funnel-note">A 使用宽选集合 A1_WIDE；主候选路径 A 使用聚焦集合 A3_FOCUSED。扩展记录独立列于附表，不作为直接互作配对。</p></section>
  <div className="filterbar">
   <div className="search-input"><Search size={16}/><Input aria-label="搜索扩展候选" placeholder="候选或种子基因、编号" value={query} onChange={e=>onQuery(e.target.value)}/></div>
   <Select value={route} onValueChange={setRoute}><SelectTrigger className="filter-select" aria-label="筛选扩展来源"><SelectValue/></SelectTrigger><SelectContent><SelectItem value="all">全部来源路径</SelectItem>{['A1','B1','B2'].map(r=><SelectItem key={r} value={r}>{r==='A1'?'A1 宽选集合':r+' 候选集合'}</SelectItem>)}</SelectContent></Select>
   <Select value={functional} onValueChange={setFunctional}><SelectTrigger className="filter-select" aria-label="筛选扩展功能证据"><SelectValue/></SelectTrigger><SelectContent><SelectItem value="all">全部功能支持</SelectItem><SelectItem value="mf">含共同分子功能</SelectItem><SelectItem value="bpmf">过程与功能联合条件</SelectItem><SelectItem value="reactome">含直接Reactome事件</SelectItem></SelectContent></Select>
   <Button variant="outline" disabled={!rows.length} onClick={()=>exportTSV('同源功能扩展.tsv',rows)}><Download size={15}/>导出{rows.length}条</Button>
  </div>
  <section className="panel homology-studio">
   <HomologyRelations rows={rows} selected={active?.candidate_accession||''} onSelect={setSelected} onProtein={onProtein}/>
   <HomologyEvidence record={active} edges={data.edges} onProtein={onProtein}/>
  </section>
  <div className="homology-table-heading"><h2>扩展候选记录</h2><Button variant="outline" aria-pressed={fullMetrics} onClick={()=>setFullMetrics(v=>!v)}><SlidersHorizontal size={15}/>{fullMetrics?'收起覆盖率':'显示覆盖率'}</Button></div>
  <DataTable data={rows} columns={cols} label="独立扩展候选"/>
  <AtlasMethods page="stage25" docs={docs}/>
 </>;
}
