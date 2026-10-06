import {createRoot} from 'react-dom/client';
import '@/app/globals.css';
import AtlasWorkspace from '@/components/atlas-workspace';

createRoot(document.getElementById('root')!).render(<AtlasWorkspace />);
