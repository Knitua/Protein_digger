import {RecordRow} from '@/lib/atlas';

export function homologyNumber(value:unknown):number|null {
 if(value==null||String(value).trim()==='')return null;
 const n=Number(value);return Number.isFinite(n)?n:null;
}
export function homologyPercent(value:unknown):string {
 const n=homologyNumber(value);return n===null?'未提供':`${(n*100).toFixed(1)}%`;
}
export function homologyTerms(value:unknown):string[] {
 return [...new Set(String(value||'').split('|').map(s=>s.trim()).filter(Boolean))];
}
/** Describe the recorded qualifying rule; never infer a measured Jaccard score. */
export function homologyFunction(record:RecordRow){
 const reason=String(record.functional_pass_reason||'');
 const bp=homologyTerms(record.shared_go_bp_direct_ids),mf=homologyTerms(record.shared_go_mf_direct_ids),events=homologyTerms(record.shared_reactome_direct_ids);
 if(reason.includes('bpmf'))return {title:'过程与功能联合支持',summary:`共同注释 ${bp.length+mf.length} 项`,rule:'GO 生物过程与分子功能联合条件：共同术语至少 3 项，IDF 加权 Jaccard ≥ 0.25。',bp,mf,events};
 if(reason.includes('go_mf'))return {title:'分子功能支持',summary:`共同分子功能 ${mf.length} 项`,rule:'GO 分子功能条件：共同术语至少 2 项，Jaccard ≥ 0.50。',bp,mf,events};
 return {title:'功能支持',summary:'查看记录的功能条件',rule:reason||'未提供判定规则',bp,mf,events};
}
