@echo off

C:
cd C:\mannual_GAMP\Tools\Batch_process_Windows
set a=%time%
echo Start time: %a%
echo start gamp_batch.py
python gamp_batch.py 2017 244 1 wum GRCE kin SF Y
set b=%time%
echo End time: %b%
pause
exit