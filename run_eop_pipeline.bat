@echo off
rem ================================================================
rem  run_eop_pipeline.bat - one-click EOP pipeline runner
rem
rem  Runs update_eop.py which will:
rem    1) download/update IERS EOP C04 (hpiers, with mirror fallback)
rem    2) merge rapid product eop\finals.daily into C04 if present
rem       (downloaded manually from https://maia.usno.navy.mil/ser7/finals.daily)
rem    3) build the C library (cmake + MSVC v145)
rem    4) run all tests (C + Python cross-validation + backtest)
rem    5) regenerate EOP_report.md
rem
rem  Usage:
rem    run_eop_pipeline.bat                 full run
rem    run_eop_pipeline.bat --skip-download skip downloads, use existing data
rem    Extra args are passed to update_eop.py
rem
rem  NOTE: keep this file ASCII/English only (GBK console safety).
rem ================================================================
setlocal EnableExtensions
cd /d "%~dp0" || exit /b 1

rem Make python output readable on the default GBK console
set "PYTHONIOENCODING=gbk"

set "RAPID_ARG="
if exist "eop\finals.daily" (
    set "RAPID_ARG=--input-rapid eop/finals.daily"
    echo [INFO] rapid product found: eop\finals.daily
) else (
    echo [INFO] eop\finals.daily not found, rapid merge disabled
    echo [INFO] manual download URL: https://maia.usno.navy.mil/ser7/finals.daily
)

where uv >nul 2>nul
if errorlevel 1 (
    echo [ERROR] uv not found in PATH. Install uv or fix PATH.
    goto :finish_fail
)

echo [INFO] running EOP pipeline ...
uv run --no-project --with numpy --with pyerfa python -I update_eop.py %RAPID_ARG% %*
set "RC=%ERRORLEVEL%"
echo.

if "%RC%"=="0" (
    echo [OK] pipeline finished successfully. Report: EOP_report.md
) else (
    echo [FAIL] pipeline exited with code %RC%. See messages above.
)
if not exist "EOP_report.md" echo [WARN] EOP_report.md was not generated.

rem Pause only when started by double-click (cmd /c), not from a terminal
echo %cmdcmdline% | find /i "/c" >nul && pause
exit /b %RC%

:finish_fail
echo %cmdcmdline% | find /i "/c" >nul && pause
exit /b 1
