preciseTarget，用于计算Cas12a2特异性识别靶标RNA底物的crRNA以及其活性和区分度

**使用方法**
1. 下载该仓库到本地并解压；
2. 使用Rstudio打开解压后的文件夹，并打开app.R文件；
3. 点击`Run App`即可启动^_^

**Tips**
若报错提示缺少依赖包，则安装所有依赖包之后再启动。


```R
# 依赖包有
depencencies <- c("shiny", "dplyr", "shinythemes", "stringr", "ggplot2", "DT", "patchwork", 
  "readxl", "RColorBrewer", "ggsci", "reshape2", "ggseqlogo")
# 安装依赖包
BiocManager::install(depencencies)
```
