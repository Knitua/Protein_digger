"use client";
import {Dna,Layers3,Waypoints,ChevronRight} from 'lucide-react';
import {RecordRow} from '@/lib/atlas';
import {homologyFunction,homologyNumber,homologyPercent} from '@/lib/homology-display';

export default function HomologyEvidence({record,edges,onProtein}:{record:RecordRow|null;edges:RecordRow[];onProtein:(id:string)=>void}){
 if(!record)return <div className="homology-comparison empty-state">没有匹配记录，请调整筛选条件。</div>;
 const evidence=homologyFunction(record),identity=homologyNumber(record.sequence_identity);
 const edge=edges.find(e=>(e.left_accession===record.best_seed_accession&&e.right_accession===record.candidate_accession)||(e.right_accession===record.best_seed_accession&&e.left_accession===record.candidate_accession));
 const weighted=String(record.functional_pass_reason).includes('bpmf');
 const functionalScore=edge?.[weighted?'go_bpmf_direct_idf_jaccard':'go_mf_direct_jaccard'];
 const metrics=[['序列一致性','sequence_identity'],['候选覆盖率','candidate_coverage'],['种子覆盖率','seed_coverage'],['结构域边界最小 IoU','domain_boundary_min_normalized_iou']] as const;
 return <div className="homology-comparison" aria-label="选中关系的支持证据">
  <header><h2>候选与种子对照</h2><span className="homology-appendix-badge">独立扩展</span></header>
  <div className="homology-pair-title">
   <button className="homology-pair-candidate" onClick={()=>onProtein(record.candidate_accession)}><span>扩展候选</span><strong>{record.candidate_gene}</strong><code>{record.candidate_accession}</code></button>
   <span className="homology-pair-line" aria-hidden="true"/>
   <button className="homology-pair-seed" onClick={()=>onProtein(record.best_seed_accession)}><span>支持种子</span><strong>{record.best_seed_gene}</strong><code>{record.best_seed_accession}</code></button>
  </div>
  <div className="homology-identity-score"><div><span><Dna size={18}/>序列一致性</span><strong>{homologyPercent(record.sequence_identity)}</strong></div><div className="homology-identity-track" aria-hidden="true"><span style={{width:identity===null?'0%':`${Math.max(0,Math.min(1,identity))*100}%`}}/></div></div>
  <div className="homology-evidence-lines">
   <div><span className="homology-evidence-icon"><Waypoints size={19}/></span><div><h3>{evidence.title}</h3><p>{[evidence.bp.length?`共同生物过程 ${evidence.bp.length} 项`:'',evidence.mf.length?`共同分子功能 ${evidence.mf.length} 项`:'',evidence.events.length?`共同直接事件 ${evidence.events.length} 项`:''].filter(Boolean).join(' · ')||'共同功能条目未提供'}</p></div></div>
   <div><span className="homology-evidence-icon"><Layers3 size={19}/></span><div><h3>{record.domain_order_copy_exact==='yes'?'结构域顺序与拷贝数一致':'结构域比较'}</h3><p>{record.domain_architecture_source||'来源未提供'} 结构域注释</p></div></div>
  </div>
  <details className="homology-full-evidence" key={record.candidate_accession}><summary><ChevronRight size={16}/>完整指标与功能条件</summary>
   <dl className="homology-detail-metrics">{metrics.map(([label,key])=><div key={key}><dt>{label}</dt><dd>{homologyNumber(record[key])===null?'未提供':record[key]}</dd></div>)}</dl>
   <h3>功能判定条件</h3><p>{evidence.rule}</p><p className="homology-detail-note">本对蛋白的{weighted?' IDF 加权 Jaccard':'分子功能 Jaccard'}：{homologyNumber(functionalScore)===null?'未提供':functionalScore}。</p>
   <dl className="homology-term-list"><dt>共同 GO 生物过程</dt><dd>{evidence.bp.join(' · ')||'未提供共同条目'}</dd><dt>共同 GO 分子功能</dt><dd>{evidence.mf.join(' · ')||'未提供共同条目'}</dd><dt>共同 Reactome 直接事件</dt><dd>{evidence.events.join(' · ')||'未提供共同条目'}</dd></dl>
  </details>
 </div>;
}
