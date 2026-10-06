"use client";
import {useEffect,useState} from 'react';
import {BookOpen,Download} from 'lucide-react';
import MethodBody from '@/components/method-body';
import {Doc} from '@/lib/release';
import {methodIds,methodSummaries} from '@/lib/atlas';
export function ReaderBody({doc}:{doc:Doc}){const [body,setBody]=useState('');useEffect(()=>{let active=true;fetch(`/Protein_digger/reader-methods/${encodeURIComponent(doc.file)}`).then(r=>{if(!r.ok)throw Error();return r.text()}).then(t=>{if(active)setBody(t)}).catch(()=>{if(active)setBody('完整方法暂未载入，请刷新后重试。')});return()=>{active=false}},[doc]);return <MethodBody text={body||'正在载入研究方法…'}/>;}
export default function AtlasMethods({page,docs}:{page:string;docs:Doc[]}){const [open,setOpen]=useState(false);return <section className="panel atlas-methods"><div className="atlas-method-heading"><BookOpen size={21}/><h2>研究方法</h2></div><p>{methodSummaries[page]}</p><details open={open} onToggle={e=>setOpen(e.currentTarget.open)}><summary>完整方法与下载</summary>{open&&docs.filter(d=>methodIds[page]?.includes(d.id)).map(d=><article key={d.id}><div className="method-file-links"><a href={`/Protein_digger/reader-methods/${encodeURIComponent(d.file)}`} download><Download size={14}/>方法 MD</a><a href={`/Protein_digger/reader-methods/${encodeURIComponent(d.file.replace(/\.md$/,'.pdf'))}`} download><Download size={14}/>方法 PDF</a></div><ReaderBody doc={d}/></article>)}</details></section>}
