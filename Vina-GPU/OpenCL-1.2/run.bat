@echo off
cd /d "%~dp0"
REM =============================================
REM  FBKdock VinaGPU — Windows Starter
REM  Cift tikla calistir
REM =============================================

REM ==== KULLANICI AYARLARI (sadece burayi degistir) ====
set LOOP=3
set METHOD=3
set CONFIG=config.txt
REM =====================================================

REM  METHOD secenekleri:
REM    1 = Tek log dosyasi (log.txt)
REM    2 = Split log (cikti klasorune)
REM    3 = Split log + logs klasoru (ONERILEN)

if not exist logs mkdir logs

for /l %%i in (1,1,%LOOP%) do (
    echo.
    echo ==========================================
    echo   Docking %%i / %LOOP%  —  %date% %time%
    echo ==========================================
    
    if %METHOD%==1 Vina-GPU.exe --config %CONFIG% --log log.txt
    if %METHOD%==2 Vina-GPU.exe --config %CONFIG% --split_log
    if %METHOD%==3 Vina-GPU.exe --config %CONFIG% --split_log --log_dir logs
    
    echo.
    echo ==========================================
    echo   Docking %%i tamamlandi
    echo ==========================================
)

echo.
echo Tum docking islemleri bitti. Cikmak icin bir tusa basin...
pause >nul
