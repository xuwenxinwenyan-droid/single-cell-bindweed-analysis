@echo off
chcp 65001 >nul
echo ========================================
echo 安装 Python 依赖包
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

echo [信息] 开始安装依赖包...
echo.
pip install -r requirements.txt

if errorlevel 1 (
    echo.
    echo [错误] 安装失败！
    echo 请检查网络连接或手动安装：
    echo pip install scanpy python-igraph leidenalg numpy pandas matplotlib seaborn
    pause
    exit /b 1
)

echo.
echo [成功] 所有依赖包已安装完成！
echo.
pause


