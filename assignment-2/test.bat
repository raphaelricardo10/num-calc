@echo off
setlocal EnableDelayedExpansion
for /f "tokens=1-4 delims=/ " %%i in ("%date%") do (
    set /a day=%%i-1
    set day=0!day!
    set day=!day:~-2!
    set month=%%j
    set year=%%k
)
set datestr=%month%%day%%year%
DEL /Q /F /S "*!datestr!*.bak"