"use client";
import {useId,type KeyboardEvent} from 'react';

/** Frozen implementation schematic. Geometry does not encode attention measurements. */
export default function NucleusArchitecture({selected,onSelect}:{selected:string;onSelect:(id:string)=>void}){
 const uid=useId().replace(/:/g,'');
 const activate=(e:KeyboardEvent<SVGGElement>,id:string)=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();onSelect(id)}};
 const group=(id:string,name:string)=>({role:'button' as const,tabIndex:0,'aria-label':name,'aria-pressed':selected===id,onClick:()=>onSelect(id),onKeyDown:(e:KeyboardEvent<SVGGElement>)=>activate(e,id),className:`nm-arch-module ${selected===id?'is-selected':''}`});
 const arrow=`url(#${uid}-arrow)`,teal=`url(#${uid}-teal-arrow)`;
 return <svg className="nm-architecture-svg" viewBox="0 0 1360 492" role="group" aria-label="Nucleus Specialist 双分支模型架构">
  <defs>
   <linearGradient id={`${uid}-background`} x1="0" y1="0" x2="1" y2="1"><stop stopColor="#f7faff"/><stop offset="1" stopColor="#f2f9f7"/></linearGradient>
   <linearGradient id={`${uid}-encoder`} x1="0" y1="0" x2="1" y2="1"><stop stopColor="#edf3ff"/><stop offset="1" stopColor="#dce9f6"/></linearGradient>
   <linearGradient id={`${uid}-specialist`} x1="0" y1="0" x2="1" y2="1"><stop stopColor="#f1faf7"/><stop offset="1" stopColor="#e7f4f0"/></linearGradient>
   <marker id={`${uid}-arrow`} markerWidth="7" markerHeight="7" refX="5" refY="3.5" orient="auto"><path d="M0 0 L6 3.5 L0 7" fill="none" stroke="#819bbb" strokeWidth="1.1" strokeLinejoin="round"/></marker>
   <marker id={`${uid}-teal-arrow`} markerWidth="7" markerHeight="7" refX="5" refY="3.5" orient="auto"><path d="M0 0 L6 3.5 L0 7" fill="none" stroke="#4e9b91" strokeWidth="1.1" strokeLinejoin="round"/></marker>
  </defs>
  <rect width="1360" height="492" rx="18" fill={`url(#${uid}-background)`}/>
  <path d="M-40 429 C210 310 282 489 533 446 S930 397 1440 470" stroke="#dbe8ef" fill="none"/>
  <path d="M-40 50 C160 159 330 -11 573 64 S1000 25 1400 54" stroke="#e2eaf4" fill="none"/>
  <g>
   <text x="24" y="133" className="nm-svg-eyebrow">序列输入</text>
   <rect x="24" y="152" width="150" height="172" rx="14" fill="#fff" stroke="#d5e2ef"/>
   <text x="99" y="181" textAnchor="middle" className="nm-svg-title">蛋白序列</text>
   {'MKALR'.split('').map((letter,i)=><g key={i}><rect x={40+i*24} y="200" width="22" height="27" rx="5" fill={i%2?'#e5f2ee':'#e8effa'}/><text x={51+i*24} y="218" textAnchor="middle" fontSize="12" fill="#486b87">{letter}</text></g>)}
   <text x="99" y="251" textAnchor="middle" className="nm-svg-copy">保留 N / C 两端</text>
   <path d="M43 271 H83 M115 271 H155" stroke="#7fadd2" strokeWidth="5" strokeLinecap="round"/><text x="99" y="276" textAnchor="middle" fill="#8fa7b9">···</text>
   <text x="99" y="306" textAnchor="middle" className="nm-svg-small">最长 1,022 残基</text>
  </g>
  <path d="M175 239 H207" className="nm-svg-link" markerEnd={arrow}/>
  <g {...group('encoder','序列表征模块')}>
   <rect className="nm-module-outline" x="216" y="152" width="188" height="172" rx="14" fill={`url(#${uid}-encoder)`} stroke="#bfcfe7"/>
   <text x="310" y="181" textAnchor="middle" className="nm-svg-title">ESM2-3B</text>
   {[0,1,2].map(i=><g key={i} transform={`translate(${257-i*8} ${200+i*10})`}><rect width="111" height="41" rx="5" fill={['#cedef2','#dce7f6','#edf3fc'][i]} stroke="#aabfdc"/>{[0,1,2,3,4,5,6].map(j=><path key={j} d={`M${12+j*14} 9 V32`} stroke="#91aecf" strokeWidth="3" opacity=".5"/>)}</g>)}
   <text x="310" y="286" textAnchor="middle" className="nm-svg-copy">残基特征 · L × 2,560</text>
   <text x="310" y="307" textAnchor="middle" className="nm-svg-small">第 36 层 · 冻结参数</text>
  </g>
  <path d="M404 239 H423 Q432 239 432 230 V111 Q432 102 441 102" className="nm-svg-link" markerEnd={arrow}/>
  <path d="M404 239 H423 Q432 239 432 248 V337 Q432 346 441 346" className="nm-svg-link is-teal" markerEnd={teal}/>
  <circle cx="422" cy="239" r="3" fill="#799ab8"/>
  <g {...group('general','通用定位分支模块')}>
   <rect className="nm-module-outline" x="444" y="22" width="596" height="190" rx="16" fill="#f6f9ff" stroke="#ccdaed"/>
   <text x="464" y="50" className="nm-svg-title">多标签定位分支</text><text x="1019" y="50" textAnchor="end" className="nm-svg-small">11 个标签查询</text>
   <rect x="461" y="77" width="115" height="53" rx="9" fill="#fff" stroke="#d1dff0"/><text x="518" y="99" textAnchor="middle" className="nm-svg-copy">特征投影</text><text x="518" y="117" textAnchor="middle" className="nm-svg-small">2,560 → 512</text>
   <rect x="461" y="145" width="115" height="46" rx="9" fill="#e9f0fb" stroke="#cfdcf0"/>
   {[0,1,2,3,4,5,6,7,8,9,10].map(i=><rect key={i} x={473+i*8} y="153" width="5" height="9" rx="1.5" fill="#7999ca"/>)}
   <text x="518" y="180" textAnchor="middle" className="nm-svg-small">可学习标签查询</text>
   <path d="M576 103 H610 M576 168 H595 Q604 168 604 159 V147 H614" className="nm-svg-link" markerEnd={arrow}/>
   <rect x="619" y="82" width="166" height="100" rx="10" fill="#fff" stroke="#bdcfe7"/>
   <text x="702" y="106" textAnchor="middle" className="nm-svg-copy">多头查询注意力</text>
   {[0,1,2,3].map(i=><g key={i}><rect x={637+i*34} y="117" width="27" height="24" rx="4" fill="#e9f0fc" stroke="#b6c9e7"/><text x={650+i*34} y="133" textAnchor="middle" fontSize="10" fill="#6683ad">{i+1}</text></g>)}
   <text x="702" y="165" textAnchor="middle" className="nm-svg-small">掩码 → 汇聚</text>
   <path d="M785 133 H813" className="nm-svg-link" markerEnd={arrow}/>
   <rect x="821" y="107" width="82" height="52" rx="9" fill="#fff" stroke="#cedcef"/><text x="862" y="129" textAnchor="middle" className="nm-svg-copy">分类头</text><text x="862" y="147" textAnchor="middle" className="nm-svg-small">逐标签</text>
   <path d="M903 133 H925" className="nm-svg-link" markerEnd={arrow}/>
   <rect x="932" y="96" width="86" height="78" rx="9" fill="#e9f0fb" stroke="#baccdf"/>
   <text x="975" y="118" textAnchor="middle" className="nm-svg-copy">通用输出</text>
   {[0,1,2,3,4,5,6,7,8,9,10].map(i=><rect key={i} x={946+i*5.5} y="131" width="3.5" height="20" rx="1" fill={i===0?'#379a91':'#91aad0'}/>)}
   <text x="975" y="166" textAnchor="middle" className="nm-svg-small">11 logits</text>
  </g>
  <g {...group('specialist','核定位专用分支模块')}>
   <rect className="nm-module-outline" x="444" y="251" width="596" height="193" rx="16" fill={`url(#${uid}-specialist)`} stroke="#bfded4"/>
   <text x="464" y="279" className="nm-svg-title">核定位专用分支</text><text x="1019" y="279" textAnchor="end" className="nm-svg-small">3 个信号约束查询</text>
   <rect x="461" y="318" width="115" height="57" rx="9" fill="#fff" stroke="#c4ddd5"/><text x="518" y="341" textAnchor="middle" className="nm-svg-copy">特征投影</text><text x="518" y="361" textAnchor="middle" className="nm-svg-small">2,560 → 384</text>
   <path d="M576 347 H582 V310 Q582 301 591 301 H800 Q809 301 809 310 V324" className="nm-svg-link is-teal" markerEnd={teal}/>
   {['NLS','NES','非核信号'].map((name,i)=><g key={name}><rect x="595" y={313+i*34} width="110" height="27" rx="6" fill="#fff" stroke="#b7d7cd"/><circle cx="609" cy={326.5+i*34} r="3" fill="#4a9f90"/><text x="657" y={331+i*34} textAnchor="middle" className="nm-svg-copy">{name}</text><path d={`M705 ${326.5+i*34} C727 ${326.5+i*34} 723 353 744 353`} className="nm-svg-link is-teal"/></g>)}
   <rect x="744" y="329" width="128" height="61" rx="10" fill="#fff" stroke="#b2d4cb"/><text x="808" y="351" textAnchor="middle" className="nm-svg-copy">四头注意力</text><text x="808" y="373" textAnchor="middle" className="nm-svg-small">掩码汇聚 · 3 组</text>
   <path d="M872 355 H891" className="nm-svg-link is-teal" markerEnd={teal}/>
   <rect x="897" y="310" width="121" height="96" rx="10" fill="#fff" stroke="#b3d6cb"/><text x="957" y="333" textAnchor="middle" className="nm-svg-copy">共享专家映射</text>
   {['NLS','NES','非核'].map((n,i)=><g key={n}><rect x={907+i*34} y="346" width="30" height="24" rx="4" fill="#e8f4ef"/><text x={922+i*34} y="361" textAnchor="middle" fontSize="9" fill="#3e7b70">{n}</text><text x={922+i*34} y="391" textAnchor="middle" className="nm-svg-copy">{i===2?'−':'+'}</text></g>)}
   <text x="809" y="429" textAnchor="middle" className="nm-svg-formula">z专用 = eNLS + eNES − e非核 + b</text>
  </g>
  <path d="M1018 133 H1055 Q1064 133 1064 142 V218 Q1064 227 1073 227 H1084" className="nm-svg-link" markerEnd={arrow}/>
  <path d="M1018 356 H1055 Q1064 356 1064 347 V294 Q1064 285 1073 285 H1084" className="nm-svg-link is-teal" markerEnd={teal}/>
  <path d="M1018 121 H1039 Q1049 121 1049 111 V89 Q1049 79 1059 79 H1217 Q1227 79 1227 89 V111 H1245" className="nm-svg-link" markerEnd={arrow}/>
  <g {...group('fusion','门控融合模块')}>
   <rect className="nm-module-outline" x="1093" y="183" width="125" height="148" rx="15" fill="#fff" stroke="#acccc9"/>
   <text x="1155" y="208" textAnchor="middle" className="nm-svg-title">门控融合</text>
   <circle cx="1155" cy="242" r="20" fill="#e8f4f1" stroke="#8cbeb3"/><text x="1155" y="249" textAnchor="middle" fontSize="22" fill="#377e73">+</text>
   <text x="1155" y="284" textAnchor="middle" className="nm-svg-copy">z通用核 + αz专用</text><text x="1155" y="312" textAnchor="middle" className="nm-svg-small">α = sigmoid(g)</text>
  </g>
  <path d="M1218 256 H1245" className="nm-svg-link is-teal" markerEnd={teal}/>
  <rect x="1253" y="221" width="84" height="74" rx="12" fill="#2f877e"/>
  <text x="1295" y="245" textAnchor="middle" fill="#fff" fontSize="14" fontWeight="600">核定位</text><text x="1295" y="268" textAnchor="middle" fill="#e5f7f2" fontSize="11">sigmoid</text><text x="1295" y="284" textAnchor="middle" fill="#e5f7f2" fontSize="10">定位概率</text>
  <rect x="1253" y="84" width="84" height="58" rx="10" fill="#edf3fb" stroke="#cfdef0"/><text x="1295" y="109" textAnchor="middle" className="nm-svg-copy">其余标签</text><text x="1295" y="130" textAnchor="middle" className="nm-svg-small">保持通用输出</text>
  <text x="310" y="465" textAnchor="middle" className="nm-svg-small">共享残基表征</text>
  <text x="741" y="470" textAnchor="middle" className="nm-svg-small">联合训练约束：分类 · 查询多样性 · 注意力平滑 · 信号约束 · 概率一致性</text>
 </svg>;
}
