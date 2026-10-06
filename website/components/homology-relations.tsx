"use client";
import {useMemo} from 'react';
import {RecordRow} from '@/lib/atlas';
import {homologyPercent} from '@/lib/homology-display';

/** Each connection is one frozen candidate–seed support record, not a PPI edge. */
export default function HomologyRelations({rows,selected,onSelect,onProtein}:{rows:RecordRow[];selected:string;onSelect:(r:RecordRow)=>void;onProtein:(id:string)=>void}){
 const groups=useMemo(()=>{
  const map=new Map<string,RecordRow[]>();
  for(const r of rows)map.set(r.best_seed_accession,[...(map.get(r.best_seed_accession)||[]),r]);
  return [...map.values()].sort((a,b)=>a[0].best_seed_gene.localeCompare(b[0].best_seed_gene)).map(rs=>rs.sort((a,b)=>a.candidate_gene.localeCompare(b.candidate_gene)||a.candidate_accession.localeCompare(b.candidate_accession)));
 },[rows]);
 return <div className="homology-map">
  <header className="homology-map-heading"><h2>同源支持关系</h2><span>{groups.length} 个支持种子 · {rows.length} 个扩展候选</span></header>
  <div className="homology-column-head" aria-hidden="true"><span>支持种子</span><span>扩展候选</span><span>序列一致性</span></div>
  {!rows.length?<div className="empty-state">当前条件下没有扩展记录</div>:<ul className="homology-groups" aria-label="按支持种子分组的同源功能关系">
   {groups.map(rs=>{const seed=rs[0];return <li key={seed.best_seed_accession} className={`homology-group ${rs.some(r=>r.candidate_accession===selected)?'has-selection':''}`}>
    <button className="homology-seed" onClick={()=>onProtein(seed.best_seed_accession)} aria-label={`查看支持种子 ${seed.best_seed_gene} ${seed.best_seed_accession} 的蛋白详情`}><strong>{seed.best_seed_gene}</strong><code>{seed.best_seed_accession}</code></button>
    <svg className="homology-branches" viewBox={`0 0 60 ${rs.length*46}`} preserveAspectRatio="none" aria-hidden="true">
     {rs.map((r,i)=><path key={r.candidate_accession} d={`M 0 ${rs.length*23} C 30 ${rs.length*23}, 30 ${i*46+23}, 60 ${i*46+23}`} className={r.candidate_accession===selected?'is-selected':''}/>)}
     <circle cx="2" cy={rs.length*23} r="3"/>
    </svg>
    <div className="homology-candidates">{rs.map(r=><button key={r.candidate_accession} className={`homology-candidate ${r.candidate_accession===selected?'is-selected':''}`} aria-pressed={r.candidate_accession===selected} aria-label={`查看 ${r.candidate_gene} 与 ${r.best_seed_gene} 的同源支持`} onClick={()=>onSelect(r)}><span><strong>{r.candidate_gene}</strong><code>{r.candidate_accession}</code></span><b>{homologyPercent(r.sequence_identity)}</b></button>)}</div>
   </li>})}
  </ul>}
  <p className="homology-map-caption">连线表示同源与功能支持，不表示直接蛋白互作。</p>
 </div>;
}
