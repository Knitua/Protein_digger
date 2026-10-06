import {defineConfig} from 'vite';
import react from '@vitejs/plugin-react';
import {fileURLToPath, URL} from 'node:url';

export default defineConfig({
  base: '/Protein_digger/',
  plugins: [react()],
  resolve: {alias: {'@': fileURLToPath(new URL('.', import.meta.url))}},
  server: {host: '127.0.0.1'},
  preview: {host: '127.0.0.1'},
  build: {sourcemap: false, chunkSizeWarningLimit: 2000},
});
