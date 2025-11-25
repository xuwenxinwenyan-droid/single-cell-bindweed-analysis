@echo off
chcp 65001 >nul
echo ========================================
echo Python 快速安装指南
echo ========================================
echo.
echo 检测到系统上的 Python 环境有问题。
echo 建议安装一个独立的 Python 环境。
echo.
echo 方法1：从官网安装（推荐）
echo   1. 访问: https://www.python.org/downloads/
echo   2. 下载 Python 3.9 或 3.10
echo   3. 安装时勾选 "Add Python to PATH"
echo   4. 安装完成后重新运行 "运行分析.bat"
echo.
echo 方法2：使用 Microsoft Store 安装
echo   1. 打开 Microsoft Store
echo   2. 搜索 "Python 3.10" 或 "Python 3.9"
echo   3. 点击安装
echo   4. 安装完成后重新运行 "运行分析.bat"
echo.
echo 方法3：使用 Anaconda（推荐用于科学计算）
echo   1. 访问: https://www.anaconda.com/products/individual
echo   2. 下载并安装 Anaconda
echo   3. 打开 Anaconda Prompt
echo   4. 运行: conda create -n scRNA python=3.9
echo   5. 运行: conda activate scRNA
echo   6. 运行: pip install -r requirements.txt
echo   7. 运行: python run_analysis.py
echo.
pause


