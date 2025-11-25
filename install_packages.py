#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
直接安装依赖包的脚本（绕过损坏的 pip）
"""
import subprocess
import sys
import os

# 需要安装的包
packages = [
    'chardet',
    'charset-normalizer',
    'numpy',
    'pandas',
    'matplotlib',
    'seaborn',
    'scanpy',
    'python-igraph',
    'leidenalg'
]

print("="*70)
print("安装 Python 依赖包")
print("="*70)
print()

# 尝试使用 pip 安装
python_exe = sys.executable

for package in packages:
    print(f"Installing {package}...")
    try:
        # 尝试使用 pip
        result = subprocess.run(
            [python_exe, '-m', 'pip', 'install', package, '--user', '--quiet'],
            capture_output=True,
            text=True,
            timeout=300
        )
        if result.returncode == 0:
            print(f"  [OK] {package} installed successfully")
        else:
            print(f"  [WARN] {package} installation failed: {result.stderr[:100]}")
    except Exception as e:
        print(f"  [ERROR] {package} installation error: {e}")

print()
print("="*70)
print("Installation complete!")
print("="*70)
print()
print("Please check the output above to confirm all packages are installed.")
print("If some packages failed, try manually:")
print("  python -m pip install <package_name> --user")

