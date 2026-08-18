#!/usr/bin/env bash
cd "$(dirname "$0")"
# =============================================
# FBKdock VinaGPU - Cross-Platform Starter
# Windows: Git Bash / WSL ile calisir
# Linux:   Native calisir
# =============================================

###### KULLANICI AYARLARI (sadece burayi degistir) ######
LOOP=3              # Kac kez pespese docking yapilsin
METHOD=3            # 1=tek log, 2=split log, 3=split+logs klasoru
CONFIG="config.txt" # Config dosyasi
##########################################################

case "$(uname -s)" in
    Linux*)  EXE="./Vina-GPU" ;;
    *)       EXE="./Vina-GPU.exe" ;;
esac

chmod +x "$EXE" 2>/dev/null || true

mkdir -p logs

for ((i=1; i<=LOOP; i++)); do
    echo ""
    echo "=========================================="
    echo "  Docking #$i / $LOOP  -  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "=========================================="
    
    case $METHOD in
        1) $EXE --config "$CONFIG" --log log.txt ;;
        2) $EXE --config "$CONFIG" --split_log ;;
        3) $EXE --config "$CONFIG" --split_log --log_dir logs ;;
    esac
    
    echo ""
    echo "=========================================="
    echo "  Docking #$i tamamlandi"
    echo "=========================================="
done
