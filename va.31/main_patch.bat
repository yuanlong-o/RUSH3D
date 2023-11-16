@echo off

REM enable local change
setlocal EnableDelayedExpansion

set /p codefile=Dir of code file:
set /p patchnum=Which patches do you want to batch process:
set /p waittime=Set delay time as:

@echo starting batch process...
@echo %patchnum%
REM run patch open matlab
for %%i in %patchnum% do (

	matlab -r "cd %codefile%/main_patch%%i, run main_patch_%%i.m"
	choice /t %waittime% /d y /n >nul
) 
pause






