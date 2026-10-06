"use client";
import {useEffect,useState} from "react";
export type RecordRow=Record<string,any>;
const requests=new Map<string,Promise<any>>();
export function readAtlas<T>(file:string):Promise<T>{const url=`/Protein_digger/atlas/${file}.json`;if(!requests.has(url))requests.set(url,fetch(url).then(r=>{if(!r.ok)throw Error('数据加载失败，请刷新后重试');return r.json()}).catch(e=>{requests.delete(url);throw e}));return requests.get(url)!;}
export function useAtlas<T>(file:string|null){const [data,setData]=useState<T|null>(null),[error,setError]=useState('');useEffect(()=>{let active=true;setData(null);setError('');if(file)readAtlas<T>(file).then(v=>{if(active)setData(v)}).catch(e=>{if(active)setError(String(e.message))});return()=>{active=false}},[file]);return {data,error};}
export const display=(v:unknown,digits=3)=>v==null||v===''?'未提供':typeof v==='number'?Number.isInteger(v)?String(v):v.toFixed(digits):String(v);
export const regulatoryLabels:Record<string,string>={KNOWN_CURATED:'人工整理调控蛋白',KNOWN_FUNCTIONAL_DIRECT:'直接调控功能证据',KNOWN_COMPLEX_ASSOCIATED:'调控复合体或功能关联',NO_KNOWN_REGULATORY_EVIDENCE:'无已知调控证据'};
export const tierLabels:Record<string,string>={KNOWN_REGULATORY:'已知调控蛋白',KNOWN_ASSOCIATED:'已知调控关联',UNDER_CHARACTERIZED:'研究不足',MIXED_MOONLIGHTING:'混合功能／潜在多功能',HK_PROBABLE:'一般核机器倾向',HK_STRONG:'一般核机器',EVIDENCE_INCOMPLETE:'证据不完整',EVIDENCE_CONFLICT:'证据冲突'};
export const depthLabels:Record<string,string>={LOW:'低',MODERATE:'中',HIGH:'高'};
export const methodIds:Record<string,string[]>={reference:['stage0'],a:['a-level1','a-level2'],b1:['b1'],b2:['b2'],b3:['b3'],pairs:['stage2'],stage25:['stage25']};
export const methodSummaries:Record<string,string>={reference:'审校人类代表序列按适用于该序列的核定位注释分区，再按调控锚点的蛋白编号和独立基因符号排除。',a:'依据人工整理名单、直接功能证据及复合体关联划分四类，再结合功能论文、证据模态和功能过程集中度进行分层。',b1:'DeepLoc与NLSExplorer分别提供核定位支持；两项分数均严格大于0.5进入候选集合。最高分窗口不等同于经过验证的核定位信号。',b2:'排除B1候选后，将注释图距离、共享种子IDF、有向网络距离和可达种子数转为百分位，等权聚合并保留前2%。',b3:'在异构体级别评价核定位，按保守分支排除代表序列与锚点重叠，再应用PLM-interact与PPLM-PPI联合互作条件。',pairs:'代表序列路线整合HI-union和官方RF2筛选集合；异构体路线采用PLM-interact与PPLM-PPI联合预测。不同证据不合成为统一强度。',stage25:'序列同源性、结构域架构与严格功能证据共同支持候选—种子关联。此类扩展不继承种子的直接互作结论。'};
