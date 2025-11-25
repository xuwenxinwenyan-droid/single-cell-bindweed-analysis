# 检查 Python 环境脚本
Write-Host "正在检查 Python 环境..." -ForegroundColor Yellow
Write-Host ""

$pythonFound = $false
$pythonPaths = @()

# 检查常见路径
$commonPaths = @(
    "python",
    "python3",
    "py",
    "$env:USERPROFILE\AppData\Local\Programs\Python\Python*\python.exe",
    "$env:USERPROFILE\anaconda3\python.exe",
    "$env:USERPROFILE\miniconda3\python.exe",
    "C:\Program Files\Python*\python.exe",
    "C:\Python*\python.exe"
)

foreach ($path in $commonPaths) {
    try {
        if ($path -match "python|py") {
            $result = & $path --version 2>&1
            if ($LASTEXITCODE -eq 0 -or $result -match "Python") {
                Write-Host "✓ 找到 Python: $path" -ForegroundColor Green
                Write-Host "  版本: $result" -ForegroundColor Green
                $pythonFound = $true
                $pythonPaths += $path
            }
        }
    } catch {
        # 忽略错误，继续查找
    }
}

# 检查环境变量
$envPaths = $env:PATH -split ';'
foreach ($envPath in $envPaths) {
    if ($envPath -match 'python|anaconda|conda') {
        $pythonExe = Join-Path $envPath "python.exe"
        if (Test-Path $pythonExe) {
            Write-Host "✓ 在 PATH 中找到: $pythonExe" -ForegroundColor Green
            $pythonFound = $true
            $pythonPaths += $pythonExe
        }
    }
}

if (-not $pythonFound) {
    Write-Host "❌ 未找到 Python 安装！" -ForegroundColor Red
    Write-Host ""
    Write-Host "请安装 Python：" -ForegroundColor Yellow
    Write-Host "1. 访问 https://www.python.org/downloads/" -ForegroundColor Cyan
    Write-Host "2. 下载 Python 3.9 或 3.10" -ForegroundColor Cyan
    Write-Host "3. 安装时勾选 'Add Python to PATH'" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "或者使用 Anaconda：" -ForegroundColor Yellow
    Write-Host "1. 访问 https://www.anaconda.com/products/distribution" -ForegroundColor Cyan
    Write-Host "2. 下载并安装 Anaconda" -ForegroundColor Cyan
} else {
    Write-Host ""
    Write-Host "找到以下 Python 安装：" -ForegroundColor Green
    foreach ($p in $pythonPaths) {
        Write-Host "  - $p" -ForegroundColor Cyan
    }
}

Write-Host ""
Write-Host "按任意键退出..."
$null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")

