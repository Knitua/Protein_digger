import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import {execFileSync} from 'node:child_process';
import os from 'node:os';
const root=path.resolve(import.meta.dirname,'..');
const m=JSON.parse(await fs.readFile(path.join(root,'data-package/manifest.json')));
const chunks=[];
for(const file of m.parts){
 const bytes=await fs.readFile(path.join(root,'data-package',file.name));
 if(crypto.createHash('sha256').update(bytes).digest('hex')!==file.sha256)throw Error(`Data-package checksum mismatch: ${file.name}`);
 chunks.push(bytes);
}
const archive=Buffer.concat(chunks);
if(crypto.createHash('sha256').update(archive).digest('hex')!==m.archive_sha256)throw Error('Archive checksum mismatch');
const scratch=await fs.mkdtemp(path.join(os.tmpdir(),'proteindigger-data-'));
const tarPath=path.join(scratch,'public.tar.gz');
try {
 await fs.writeFile(tarPath,archive);
 const listing=execFileSync('tar',['-tzf',tarPath],{encoding:'utf8'});
 if(listing.split('\n').some(p=>p.startsWith('/')||p.split('/').includes('..')))throw Error('Unsafe archive path');
 const output=path.join(root,'public');await fs.mkdir(output,{recursive:true});
 execFileSync('tar',['-xzf',tarPath,'-C',output]);
 console.log('Frozen public data restored and archive checksum verified.');
} finally {
 await fs.unlink(tarPath);await fs.rmdir(scratch);
}
