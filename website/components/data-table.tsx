"use client";
import {useEffect,useState} from "react";
import {ColumnDef,SortingState,flexRender,getCoreRowModel,getPaginationRowModel,getSortedRowModel,useReactTable} from "@tanstack/react-table";
import {Table,TableBody,TableCell,TableHead,TableHeader,TableRow} from "@/components/ui/table";
import {Button} from "@/components/ui/button";
import {ChevronLeft,ChevronRight,ArrowUpDown} from "lucide-react";
import {fmt} from "@/lib/release";
export default function DataTable<T>({data,columns,label}:{data:T[];columns:ColumnDef<T>[];label:string}){
 const [sorting,setSorting]=useState<SortingState>([]);
 const table=useReactTable({data,columns,state:{sorting},onSortingChange:setSorting,getCoreRowModel:getCoreRowModel(),getSortedRowModel:getSortedRowModel(),getPaginationRowModel:getPaginationRowModel(),initialState:{pagination:{pageSize:25}}});
 useEffect(()=>{table.setPageIndex(0)},[data,table]);
 return <div className="table-surface"><Table aria-label={label}><TableHeader>{table.getHeaderGroups().map(g=><TableRow key={g.id}>{g.headers.map(h=><TableHead key={h.id}>{h.isPlaceholder?null:<button className="sort-head" disabled={!h.column.getCanSort()} onClick={h.column.getToggleSortingHandler()}>{flexRender(h.column.columnDef.header,h.getContext())}{h.column.getCanSort()&&<ArrowUpDown size={12}/>}</button>}</TableHead>)}</TableRow>)}</TableHeader><TableBody>{table.getRowModel().rows.length?table.getRowModel().rows.map(r=><TableRow key={r.id}>{r.getVisibleCells().map(c=><TableCell key={c.id}>{flexRender(c.column.columnDef.cell,c.getContext())}</TableCell>)}</TableRow>):<TableRow><TableCell colSpan={columns.length}><div className="empty-state">没有匹配记录<span>请调整基因、蛋白编号、筛选路径或证据类型。检索范围为当前工作区的冻结数据。</span></div></TableCell></TableRow>}</TableBody></Table><div className="pagination"><span>共 {fmt(data.length)} 条 · 每页 25 条</span><div><Button aria-label="上一页" variant="outline" size="icon-sm" disabled={!table.getCanPreviousPage()} onClick={()=>table.previousPage()}><ChevronLeft/></Button><span>{data.length?table.getState().pagination.pageIndex+1:0} / {table.getPageCount()}</span><Button aria-label="下一页" variant="outline" size="icon-sm" disabled={!table.getCanNextPage()} onClick={()=>table.nextPage()}><ChevronRight/></Button></div></div></div>;
}
