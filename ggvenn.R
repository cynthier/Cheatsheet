library("ggvenn")
ggvenn(
  data = genes,         # 数据列表, LIST
  columns = NULL,           # 对选中的列名绘图，最多选择4个，NULL为默认全选
  show_elements = F,        # 当为TRUE时，显示具体的交集情况，而不是交集个数
  label_sep = "\n",         # 当show_elements = T时生效，分隔符 \n 表示的是回车的意思
  show_percentage = T,      # 显示每一组的百分比
  digits = 1,               # 百分比的小数点位数
  fill_color = c("#E41A1C", "#1E90FF", "#FF8C00", "#80FF00"), # 填充颜色
  fill_alpha = 0.5,         # 填充透明度
  stroke_color = "white",   # 边缘颜色
  stroke_alpha = 0.5,       # 边缘透明度
  stroke_size = 0.5,        # 边缘粗细
  stroke_linetype = "solid", # 边缘线条 # 实线：solid  虚线：twodash longdash 点：dotdash dotted dashed  无：blank
  set_name_color = "black", # 组名颜色
  set_name_size = 6,        # 组名大小
  text_color = "black",     # 交集个数颜色
  text_size = 4             # 交集个数文字大小
)
