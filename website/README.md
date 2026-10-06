# ProteinDigger research workbench

**[打开完整研究工作台 / Open the workbench](https://knitua.github.io/Protein_digger/)**

[候选蛋白库](https://knitua.github.io/Protein_digger/#view=candidates) · [蛋白互作网络](https://knitua.github.io/Protein_digger/#view=pairs) · [核定位模型与评估](https://knitua.github.io/Protein_digger/#view=nucleus-model) · [方法与下载](https://knitua.github.io/Protein_digger/#view=downloads)

## 内容与边界

这是完整本地工作台的静态公开版，保留各阶段数据库、二维/三维网络、蛋白详情、搜索筛选、当前结果导出、方法 MD/PDF 与原始交付 Excel。页面不连接科研服务器，不提交预测，不收集候选查询，不包含分析遥测。GitHub Pages 负责静态文件托管；打开网页不需要运行本地服务或登录。

- 主候选集合：1,580 条路径记录，1,573 个唯一基因，1,620 个蛋白编号。
- 互作：5,957 条来源记录，5,954 个唯一蛋白对。没有提供的模型分数仍为空值。
- B3：112 个基因、152 条异构体，保留保守分支及历史 B1 排除集边界。
- 同源功能扩展：9 条独立记录，采用 v3 输入来源与支持证据；不将支持种子的直接互作结论转移给扩展候选。
- 核定位模型：独立模型研究，不替换当前 B1 模型。完整下载保留八项指标；图中展示选定指标，误差线为五折标准差。

公开版本未改变候选、配对或科研评分。只删除了未使用的旧详情副本和非发布资料的分发引用，并将两份基准元数据中的服务器路径替换为来源文件名；原始及公开副本的 SHA256 对照记录在 `public/publication-manifest.json`。科研目录、原始 Excel 和本地技术审计未改写。

原始 Excel SHA256：
`89c5d398d3bbd825f4541f34b09c5f3ecf930f0d9d2488cc15367af2f3ebc4e1`

各数据库和模型的来源见页面方法与源文献；第三方资源仍受各自许可和使用条款约束。公开浏览不意味着获得重新分发第三方原始数据库的许可。

## 本地构建

使用 Node.js 22.13 或更高版本：

```bash
cd website
npm ci
npm run prepare:data
npm test
npm run build
npm run preview
```

预览地址为 `http://127.0.0.1:4173/Protein_digger/`。数据在浏览器侧处理，大型网络与蛋白详情按需加载；首次进入完整阶段数据库可能需要下载相应数据包。网络空间布局不是分子位置，密度图也不是功能活性或富集程度。

## 发布与更新

仓库 Pages 使用 GitHub Actions。`main` 分支的 `website/` 或 Pages 工作流更新会依次恢复冻结数据、运行数据校验、TypeScript 检查和 Vite 构建，只有通过后才发布 `website/dist`。也可在 Actions → Publish research workbench → Run workflow 手动重新部署。

完整静态资源以 `data-package/` 中的分块压缩包保存，`manifest.json` 记录每块及完整归档的 SHA256。`npm run prepare:data` 校验并解包至 `public/`；并不重新生成科研结果。采用此分发形式是为了减少大批 JSON 文件的传输体积，网页仍能访问全部阶段数据、蛋白详情和下载资料。

站点子路径固定为 `/Protein_digger/`；书签采用 `#view=...`，刷新和直接访问无需服务器重写。修改科研结果需另行完成来源核验，不能通过更改显示计数绕过数据测试。发布前检查 `publication-manifest.json` 与下载文件校验值；禁止提交凭证、运行日志、服务器备份、模型权重或个人配置。
