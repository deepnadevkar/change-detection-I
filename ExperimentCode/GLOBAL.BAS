Attribute VB_Name = "Module1"
Option Explicit

Global listenToFormPecks As Boolean

Global gsProfile As String
Global gsResultsFile As String
Global gsImagePath As String
Global gbRecalculateJuicePellets As Boolean

'Global stimulus array will hold references to the
'Lead pictures on frmDisplay
Global stimulus_array(1 To 16) As Variant


'screen size
Global Const SCREEN_TOP = 0
Global Const SCREEN_LEFT = 0
Global Const SCREEN_HEIGHT = 480 '600 '480
Global Const SCREEN_WIDTH = 640 '800 '640

'LEFT    TOP -- got these from excel conversions; adjusted below manually
'365 193
'315 279.6025404
'215 279.6025404
'165 193
'215 106.3974596
'315 106.3974596







'Stimuli size and location
Global Const STIMULI_TOP = 183 '-15 '138 '198  '128 '329
Global Const STIMULI_LEFT = 385 '100 ' 300 '397
Global Const STIMULI_HEIGHT = 94 '130 = 5 cm    '90 va
Global Const STIMULI_WIDTH = 125 '166 = 6 cm    '134 va

'Stimuli size and location
Global Const LEAD1_TOP = 269.6025404 '-15 '198  '128 '329
Global Const LEAD1_LEFT = 335 '220 '397
Global Const LEAD1_HEIGHT = 94 '130 = 5 cm    '90 va
Global Const LEAD1_WIDTH = 125 '166 = 6 cm    '134 va

'Stimuli size and location
Global Const LEAD2_TOP = 269.6025404 '-15 '198  '128 '329
Global Const LEAD2_LEFT = 235 '340 '397
Global Const LEAD2_HEIGHT = 94 '130 = 5 cm    '90 va
Global Const LEAD2_WIDTH = 125 '166 = 6 cm    '134 va

'Stimuli size and location
Global Const LEAD3_TOP = 183 '-15 '198  '128 '329
Global Const LEAD3_LEFT = 185 '460 '397
Global Const LEAD3_HEIGHT = 94 '130 = 5 cm    '90 va
Global Const LEAD3_WIDTH = 125 '166 = 6 cm    '134 va

'Stimuli size and location
Global Const LEAD4_TOP = 96.3974596 '105 '198  '128 '329
Global Const LEAD4_LEFT = 235 '100 '397
Global Const LEAD4_HEIGHT = 94 '130 = 5 cm    '90 va
Global Const LEAD4_WIDTH = 125 '166 = 6 cm    '134 va

'Stimuli size and location
Global Const LEAD5_TOP = 96.3974596 '105 '198  '128 '329
Global Const LEAD5_LEFT = 335 '220 '397
Global Const LEAD5_HEIGHT = 94 '130 = 5 cm    '90 va
Global Const LEAD5_WIDTH = 125 '166 = 6 cm    '134 va

'Stimuli size and location
Global Const LEAD6_TOP = 175 '105 '198  '128 '329
Global Const LEAD6_LEFT = 250 '340 '397
Global Const LEAD6_HEIGHT = 94 '130 = 5 cm    '90 va
Global Const LEAD6_WIDTH = 125 '166 = 6 cm    '134 va

'Stimuli size and location
Global Const LEAD7_TOP = 122 '105 '198  '128 '329
Global Const LEAD7_LEFT = 525 '460 '397
Global Const LEAD7_HEIGHT = 94 '130 = 5 cm    '90 va
Global Const LEAD7_WIDTH = 125 '166 = 6 cm    '134 va

'Stimuli size and location
Global Const LEAD8_TOP = 254 '225 '198  '128 '329
Global Const LEAD8_LEFT = 75 '100 '397
Global Const LEAD8_HEIGHT = 94 '130 = 5 cm    '90 va
Global Const LEAD8_WIDTH = 125 '166 = 6 cm    '134 va

'Stimuli size and location
Global Const LEAD9_TOP = 254 '225 '198  '128 '329
Global Const LEAD9_LEFT = 225 '220 '397
Global Const LEAD9_HEIGHT = 94 '130 = 5 cm    '90 va
Global Const LEAD9_WIDTH = 125 '166 = 6 cm    '134 va

'Stimuli size and location
Global Const LEAD10_TOP = 254 '225 '198  '128 '329
Global Const LEAD10_LEFT = 375 '340 '397
Global Const LEAD10_HEIGHT = 94 '130 = 5 cm    '90 va
Global Const LEAD10_WIDTH = 125 '166 = 6 cm    '134 va

'Stimuli size and location
Global Const LEAD11_TOP = 254 '225 '198  '128 '329
Global Const LEAD11_LEFT = 525 '460 '397
Global Const LEAD11_HEIGHT = 94 '130 = 5 cm    '90 va
Global Const LEAD11_WIDTH = 125 '166 = 6 cm    '134 va

'Stimuli size and location
Global Const LEAD12_TOP = 386 '345 '198  '128 '329
Global Const LEAD12_LEFT = 75 '100 '397
Global Const LEAD12_HEIGHT = 94 '130 = 5 cm    '90 va
Global Const LEAD12_WIDTH = 125 '166 = 6 cm    '134 va

'Stimuli size and location
Global Const LEAD13_TOP = 386 '345 '198  '128 '329
Global Const LEAD13_LEFT = 225 '220 '397
Global Const LEAD13_HEIGHT = 94 '130 = 5 cm    '90 va
Global Const LEAD13_WIDTH = 125 '166 = 6 cm    '134 va

'Stimuli size and location
Global Const LEAD14_TOP = 386 '345 '198  '128 '329
Global Const LEAD14_LEFT = 375 '340 '397
Global Const LEAD14_HEIGHT = 94 '130 = 5 cm    '90 va
Global Const LEAD14_WIDTH = 125 '166 = 6 cm    '134 va

'Stimuli size and location
Global Const LEAD15_TOP = 386 '345 '198  '128 '329
Global Const LEAD15_LEFT = 525 '460 '397
Global Const LEAD15_HEIGHT = 94 '130 = 5 cm    '90 va
Global Const LEAD15_WIDTH = 125 '166 = 6 cm    '134 va

'Probe size and location
Global Const PROBE_TOP = 178 '435 '70   '450
Global Const PROBE_LEFT = 287 '400 '397  '319
Global Const PROBE_HEIGHT = 94  '90 = 3.5 cm
Global Const PROBE_WIDTH = 125 '134 = 5 cm

'Diff size and location
Global Const DIFF_TOP = 40 '70  '490
Global Const DIFF_LEFT = 143 '286 '490
Global Const DIFF_HEIGHT = 0 '56   '50 = 2cm
Global Const DIFF_WIDTH = 0 '70c    '68
'Colors
Global Const COLOR_LIGHT_RED = &H8080FF
Global Const COLOR_RED = &HC0&
Global Const COLOR_DARK_RED = &H80&
Global Const COLOR_GREEN = &HC000&
Global Const COLOR_DARK_GREEN = &H8000&
Global Const COLOR_LIGHT_GREEN = &H80FF80

Global Const COLOR_WHITE = &HFFFFFF
Global Const COLOR_BLACK = &H0&
Global Const COLOR_GRAY = &HC0C0C0

Declare Sub agOutp Lib "apigid32.dll" (ByVal portid%, ByVal outval%)
Declare Function agInp Lib "apigid32.dll" (ByVal portid%) As Integer

Declare Function ClipCursor Lib "User32" (lpRect As Any) As Long
Declare Function ClipCursorByNum Lib "User32" Alias "ClipCursor" (lpRect As Any) As Long
Declare Function GetTickCount Lib "kernel32" () As Long
Declare Function ShowCursor Lib "User32" (ByVal bShow As Integer) As Integer
Declare Function SetCursorPos Lib "User32" (ByVal X As Long, ByVal Y As Long) As Long

Type RECT
        Left As Long
        Top As Long
        Right As Long
        Bottom As Long
End Type

Public Sub gsubWait(iSeconds As Integer)

    Dim vStart
    
    vStart = Timer   ' Set start time.
    Do While Timer < vStart + iSeconds
        DoEvents    ' Yield to other processes.
    Loop
    
End Sub
Public Sub cursor_visible(mtest As Boolean)
    Dim pos As Integer
    Select Case mtest
        Case True
            Do
                pos = ShowCursor(True)
            Loop Until pos = 1
         Case False
            Do
                pos = ShowCursor(False)
            Loop Until pos = -1
    End Select
End Sub
Public Sub gsubScreenDimensions(frm As Form, iOffsetX As Integer, iOffsetY As Integer)
    
    '800 x 600 '640 by 480
    frm.Top = iOffsetY * Screen.TwipsPerPixelY
    frm.Left = iOffsetX * Screen.TwipsPerPixelX
    frm.Height = 9025 '7220
    frm.Width = 12037 '9630

End Sub


Public Sub gWaitMilliseconds(lMilliseconds As Long)
    
    Dim lStart As Long
    
    lStart = GetTickCount()   ' Set start time.
    Do While GetTickCount() < lStart + lMilliseconds
        DoEvents    ' Yield to other processes.
    Loop

End Sub

Sub Main()

    On Error GoTo Error_Main
    
    frmMain.Show
    Exit Sub
    
Error_Main:

    MsgBox Err & ":  " & Error
    Exit Sub

End Sub

Private Sub subReleaseMouse()

    Dim lpRect As RECT
    Dim lRet As Long
    
    lpRect.Top = 0
    lpRect.Left = 0
    lpRect.Right = Screen.Width / Screen.TwipsPerPixelX
    lpRect.Bottom = Screen.Height / Screen.TwipsPerPixelY
    lRet = ClipCursor(lpRect)

End Sub

Private Sub subRestrictMouse()

 Dim lRet As Long
    Dim lpRect As RECT
    
    lpRect.Top = 0
    lpRect.Left = 0
    lpRect.Right = SCREEN_WIDTH
    lpRect.Bottom = SCREEN_HEIGHT
    lRet = ClipCursor(lpRect)

End Sub

