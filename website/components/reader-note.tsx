import {type ReactNode} from 'react';
import {ChevronDown,Info} from 'lucide-react';

/** Keep scientific qualifications accessible without repeating instructions in the main view. */
export default function ReaderNote({title='证据说明',children}:{title?:string;children:ReactNode}){
 return <details className="reader-note"><summary><Info size={14}/><span>{title}</span><ChevronDown size={14}/></summary><div className="reader-note-body">{children}</div></details>;
}
