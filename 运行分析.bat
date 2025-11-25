@echo off
chcp 65001 >nul
echo ========================================
echo 单细胞RNA-seq分析脚本
echo ========================================
echo.

REM 检查 Python 是否安装
python --version >nul 2>&1
if errorlevel 1 (
    echo [错误] 未找到 Python！
    echo.
    echo 请先安装 Python：
    echo 1. 访问 https://www.python.org/downloads/
    echo 2. 下载并安装 Python 3.8 或更高版本
    echo 3. 安装时勾选 "Add Python to PATH"
    echo.
    pause
    exit /b 1
)

echo [信息] 找到 Python
python --version
echo.

REM 检查是否安装了必需的包
echo [信息] 检查必需的包...
python -c "import scanpy" >nul 2>&1
if errorlevel 1 (
    echo [警告] 缺少必需的包，正在安装...
    echo.
    pip install -r requirements.txt
    if errorlevel 1 (
        echo [错误] 包安装失败！
        echo 请手动运行: pip install -r requirements.txt
        pause
        exit /b 1
    )
    echo.
    echo [成功] 包安装完成
) else (
    echo [成功] 所有必需的包已安装
)

echo.
echo ========================================
echo 开始运行分析...
echo ========================================
echo.
echo 注意：完整分析可能需要 30-45 分钟
echo 请耐心等待...
echo.

REM 运行分析脚本
python run_analysis.py

if errorlevel 1 (
    echo.
    echo [错误] 分析过程中出现错误！
    pause
    exit /b 1
)

echo.
echo ========================================
echo 分析完成！
echo ========================================
echo.
echo 结果已保存在 ./results/ 文件夹中
echo.
pause


