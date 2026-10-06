"use client";
import {useEffect,useId,useState} from "react";
import ReactMarkdown from "react-markdown";
import remarkGfm from "remark-gfm";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css";
function Mermaid({code}:{code:string}) {
 const id=useId().replace(/[^a-zA-Z0-9]/g,"");const [svg,setSvg]=useState("");const [error,setError]=useState(false);
 useEffect(()=>{let live=true;import("mermaid").then(async({default:m})=>{
 m.initialize({startOnLoad:false,securityLevel:"strict",theme:"base",themeVariables:{primaryColor:"#eef3ff",primaryTextColor:"#18314f",primaryBorderColor:"#b7c9e8",lineColor:"#7088a8",fontFamily:"Arial, PingFang SC, sans-serif",fontSize:"14px"},flowchart:{htmlLabels:false,useMaxWidth:true}});
 try{const r=await m.render(`graph${id}`,code);if(live)setSvg(r.svg);}catch{if(live)setError(true);}
 });return()=>{live=false};},[code,id]);
 return error?<pre className="diagram-error">流程图源码（渲染失败）{code}</pre>:<div className="diagram" aria-label="方法流程图" dangerouslySetInnerHTML={{__html:svg||"正在加载流程图…"}}/>;
}
export default function MethodBody({text}:{text:string}) {
 return <article className="method-prose"><ReactMarkdown remarkPlugins={[remarkGfm,remarkMath]} rehypePlugins={[rehypeKatex]} components={{
 code({className,children,...props}){const code=String(children).replace(/\n$/,"");return className==="language-mermaid"?<Mermaid code={code}/>:<code className={className} {...props}>{children}</code>;},
 pre({children}){return <div className="code-block">{children}</div>;},
 table({children}){return <div className="method-table"><table>{children}</table></div>;},
 a({href,children}){return <a href={href} target="_blank" rel="noreferrer">{children}</a>;},
 }}>{text}</ReactMarkdown></article>;
}
