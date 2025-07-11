import matplotlib.pyplot as plt
import numpy as np

# --- 1. 准备你的数据 (用你自己的数据替换这里的示例数据) ---
# 横坐标 N (展开结束N)
N_data = np.array([0, 2, 4])  # 示例数据：从1到10的整数

# 纵坐标 F (自能能F)
# 示例数据：一些随机值，你可以替换成你的计算结果
F_data = np.array(
    [-25.49731, -25.49731 - 0.045835, -25.49731 - 0.045835 + 0.0000197534]
)
F_err = np.array([0.0, 0, 0])
# F_data = np.array([-25.49731, -25.49731 - 0.045835])

# --- 2. 创建图形 ---
# 创建一个图形和一个坐标轴对象
fig, ax = plt.subplots(figsize=(8, 6))

# --- 3. 绘制数据 ---
# 绘制折线图，并添加标记点(marker='o')以便更清晰地查看数据点
# ax.errorbar(N_data, F_data, F_err, pt="o", linestyle="--")
ax.errorbar(N_data, F_data / F_data[0], F_err / F_data[0])

# --- 4. 设置图表属性 ---
# 设置横坐标标签
ax.set_xlabel(r"N")

# 设置纵坐标标签
ax.set_ylabel(r"F")


# 显示图例 (如果有多条线，这会很有用)
ax.legend()

# --- 5. 显示图形 ---
# 自动调整布局，防止标签重叠
plt.tight_layout()

# 显示最终的图表
plt.show()
