"use client";
import {useEffect} from 'react';
import {Release} from './release';
type Registry={registerTool(tool:{name:string;title:string;description:string;inputSchema:object;annotations:object;execute:(input:unknown)=>unknown},options:{signal:AbortSignal}):void|Promise<void>};
export function useReadonlyTools(release:Release|null){
 useEffect(()=>{
  const context=(document as Document&{modelContext?:Registry}).modelContext;
  if(!release||!context?.registerTool)return;
  const lifecycle=new AbortController();
  try{void Promise.resolve(context.registerTool({
   name:'search_delivered_candidates',title:'检索冻结交付候选',
   description:'只读检索当前候选主表；不提交计算、不修改数据，不搜索同源功能扩展独立附表。最多返回 25 条及总匹配数。',
   inputSchema:{type:'object',properties:{query:{type:'string'},route:{type:'string',enum:['all','A','B1','B2','B3']}},required:['query'],additionalProperties:false},
   annotations:{readOnlyHint:true,untrustedContentHint:true},
   execute(input){
    if(!input||typeof input!=='object'||Array.isArray(input))throw Error('需要 query 对象');
    const v=input as Record<string,unknown>;
    if(typeof v.query!=='string'||Object.keys(v).some(k=>!['query','route'].includes(k))||v.route!==undefined&&!(typeof v.route==='string'&&['all','A','B1','B2','B3'].includes(v.route)))throw Error('无效查询参数');
    const q=v.query.trim().toLowerCase(),r=v.route||'all';
    const matches=release.candidates.filter(c=>(r==='all'||c.route===r)&&`${c.gene} ${c.accessions.join(' ')} ${c.name}`.toLowerCase().includes(q));
    return {release:release.release,total:matches.length,records:matches.slice(0,25)};
   }
  },{signal:lifecycle.signal})).catch(e=>console.warn('本地只读工具未注册',String(e)))}catch(e){console.warn('浏览器未支持只读工具',String(e))}
  return()=>lifecycle.abort();
 },[release]);
}
