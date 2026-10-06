export type ModelMetric={mean:number;std:number;folds:number[]};
export type BenchmarkMetric={id:string;label:string;definition:string;deeploc:ModelMetric;specialist:ModelMetric;delta_pp:number};
export type BenchmarkFold={id:number;test:number;held_in:number;inner_validation:number;train:number;thresholds:{deeploc:number;specialist:number}};
export type NucleusBenchmark={schema_version:1;name:string;baseline:string;n_samples:number;fold_count:number;ensemble:boolean;std_ddof:number;metrics:BenchmarkMetric[];folds:BenchmarkFold[];scientific_boundary:string;source_checksums:{id:string;sha256:string;download:string}[];source_commit:string;public_benchmark_commit:string;dataset_sha256:string;oof_sha256:string};
export function validateNucleusBenchmark(value:unknown):NucleusBenchmark{
 const x=value as NucleusBenchmark;
 if(x?.schema_version!==1||x.n_samples!==28300||x.fold_count!==5||x.ensemble!==false||x.std_ddof!==0||!Array.isArray(x.metrics)||x.metrics.length!==8||!Array.isArray(x.folds)||x.folds.length!==5)throw Error('模型评估数据格式不完整，请重新载入。');
 for(const m of x.metrics){for(const model of ['deeploc','specialist'] as const){const v=m[model];if(!v||!Number.isFinite(v.mean)||!Number.isFinite(v.std)||!Array.isArray(v.folds)||v.folds.length!==5||!v.folds.every(Number.isFinite))throw Error('五折评估数据不完整。')}if(!Number.isFinite(m.delta_pp))throw Error('指标差值缺失。')}
 if(!x.metrics.some(m=>m.id==='nucleus-recall'))throw Error('评估缺少召回率。');
 return x;
}
export const benchmarkBase='/Protein_digger/models/nucleus-specialist';
export const score=(n:number)=>n.toFixed(4);
export const delta=(n:number)=>`${n>0?'+':n<0?'−':''}${Math.abs(n).toFixed(2)}`;
