"use client";
import {useEffect,useState} from 'react';
const key='ep-background-motion';
export function useBackgroundMotion(){
 const [motion,setState]=useState(false);
 useEffect(()=>{const media=matchMedia('(prefers-reduced-motion: reduce)');const sync=()=>setState(!media.matches&&localStorage.getItem(key)!=='paused');sync();media.addEventListener('change',sync);window.addEventListener(key,sync);return()=>{media.removeEventListener('change',sync);window.removeEventListener(key,sync)}},[]);
 const setMotion=(value:boolean|((previous:boolean)=>boolean))=>{const next=typeof value==='function'?value(motion):value;localStorage.setItem(key,next?'enabled':'paused');window.dispatchEvent(new Event(key))};
 return [motion,setMotion] as const;
}
