import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
const root=path.resolve(import.meta.dirname,'../public');
const read=p=>fs.readFile(path.join(root,p));
const json=async p=>JSON.parse(await read(p));
const sha=b=>crypto.createHash('sha256').update(b).digest('hex');
const release=await json('data/release.json');
const {candidates,pairs}=release;
assert.equal(candidates.length,1580);
assert.equal(new Set(candidates.flatMap(c=>c.gene.split(/[;,]/).map(x=>x.trim()).filter(Boolean))).size,1573);
assert.equal(new Set(candidates.flatMap(c=>c.accessions)).size,1620);
assert.equal(pairs.length,5957);
assert.equal(new Set(pairs.map(p=>[p.candidate,p.anchor].sort().join('|'))).size,5954);
assert.equal(release.appendix.length,9);
for(const [route,count]of Object.entries({A:954,B1:429,B2:85,B3:112}))assert.equal(candidates.filter(c=>c.route===route).length,count);
assert.equal(new Set(candidates.filter(c=>c.route==='B3').flatMap(c=>c.accessions)).size,152);
assert.equal(sha(await read('downloads/第一批候选基因与对应蛋白.xlsx')),'89c5d398d3bbd825f4541f34b09c5f3ecf930f0d9d2488cc15367af2f3ebc4e1');
for(const [file,count]of Object.entries({reference:20416,annotation:3418,prediction:14625,network:14625,isoforms:22131}))assert.equal((await json(`atlas/${file}.json`)).length,count,file);
assert.equal((await json('atlas/extension.json')).records.length,9);
for(const p of pairs){assert(candidates.some(c=>c.route===p.route&&c.accessions.includes(p.candidate)),`${p.id}: candidate`);if(p.route!=='B3'){assert.equal(p.plm,null);assert.equal(p.pplm,null);}}
const benchmark=await json('models/nucleus-specialist/benchmark-v1.json');
assert.equal(benchmark.metrics.length,8);assert.equal(benchmark.n_samples,28300);
for(const m of benchmark.metrics)for(const model of ['deeploc','specialist']){const s=m[model];assert.equal(s.folds.length,5);const mean=s.folds.reduce((a,b)=>a+b,0)/5;assert(Math.abs(mean-s.mean)<1e-10);const sd=Math.sqrt(s.folds.reduce((a,b)=>a+(b-mean)**2,0)/5);assert(Math.abs(sd-s.std)<1e-10);}
for(const s of benchmark.source_checksums)assert.equal(sha(await read('models/nucleus-specialist/'+s.download)),s.sha256);
const manifest=await json('publication-manifest.json');
for(const f of manifest.files){const b=await read(f.file);assert.equal(b.length,f.bytes,f.file);assert.equal(sha(b),f.sha256,f.file);assert(b.length<100*1024*1024,'GitHub single-file limit');if(/\.(json|md|tsv)$/.test(f.file))assert(!/\/root\/|\/Users\/|autodl-tmp/.test(b.toString()),`${f.file}: internal path`);}
for(const d of await json('atlas/documents.json')){await read('reader-methods/'+d.file);await read('reader-methods/'+d.file.replace(/\.md$/,'.pdf'));}
console.log(`Verified ${manifest.files.length} published files; 1,580 candidate rows, 5,957 source rows, 5,954 unique pairs; Excel SHA256 unchanged.`);
