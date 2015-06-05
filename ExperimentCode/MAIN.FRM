VERSION 5.00
Object = "{F9043C88-F6F2-101A-A3C9-08002B2F49FB}#1.2#0"; "COMDLG32.OCX"
Begin VB.Form frmMain 
   BackColor       =   &H00808000&
   BorderStyle     =   1  'Fixed Single
   Caption         =   "Pigeon Change Detection"
   ClientHeight    =   7335
   ClientLeft      =   1875
   ClientTop       =   2250
   ClientWidth     =   10950
   ControlBox      =   0   'False
   LinkTopic       =   "Form1"
   MaxButton       =   0   'False
   MinButton       =   0   'False
   PaletteMode     =   1  'UseZOrder
   ScaleHeight     =   7335
   ScaleWidth      =   10950
   Begin MSComDlg.CommonDialog diaResults 
      Left            =   1800
      Top             =   750
      _ExtentX        =   847
      _ExtentY        =   847
      _Version        =   393216
   End
   Begin MSComDlg.CommonDialog diaFile 
      Left            =   1800
      Top             =   180
      _ExtentX        =   847
      _ExtentY        =   847
      _Version        =   393216
   End
   Begin VB.Label Label2 
      Alignment       =   2  'Center
      BackColor       =   &H00808000&
      Caption         =   "Juice Support, New Labels"
      BeginProperty Font 
         Name            =   "Times New Roman"
         Size            =   48
         Charset         =   0
         Weight          =   700
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Height          =   2295
      Left            =   360
      TabIndex        =   1
      Top             =   4920
      Width           =   10455
   End
   Begin VB.Label Label1 
      Alignment       =   2  'Center
      BackColor       =   &H00808000&
      Caption         =   "Change Detection DRO (PCI) Automatic VPD"
      BeginProperty Font 
         Name            =   "Times New Roman"
         Size            =   48
         Charset         =   0
         Weight          =   700
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      ForeColor       =   &H00000080&
      Height          =   3495
      Left            =   1320
      TabIndex        =   0
      Top             =   1560
      Width           =   9375
   End
   Begin VB.Menu mnuFile 
      Caption         =   "&File"
      Begin VB.Menu mnuNew 
         Caption         =   "&New Session"
      End
      Begin VB.Menu mnuOpen 
         Caption         =   "&Open Session"
      End
      Begin VB.Menu mnuClose 
         Caption         =   "&Close Session"
         Enabled         =   0   'False
      End
      Begin VB.Menu mnuSep2 
         Caption         =   "-"
      End
      Begin VB.Menu mnuNewResults 
         Caption         =   "N&ew Results"
      End
      Begin VB.Menu mnuOpenResults 
         Caption         =   "O&pen Results"
      End
      Begin VB.Menu mnuCloseResults 
         Caption         =   "C&lose Results"
         Enabled         =   0   'False
      End
      Begin VB.Menu mnuSep1 
         Caption         =   "-"
      End
      Begin VB.Menu mnuExit 
         Caption         =   "E&xit"
      End
   End
   Begin VB.Menu mnuDefaults 
      Caption         =   "&Defaults"
      Enabled         =   0   'False
   End
   Begin VB.Menu mnuSetup 
      Caption         =   "&Setup"
      Enabled         =   0   'False
   End
   Begin VB.Menu mnuRun 
      Caption         =   "&Run"
      Enabled         =   0   'False
   End
   Begin VB.Menu mnuResults 
      Caption         =   "&Results"
      Enabled         =   0   'False
   End
   Begin VB.Menu mnuView 
      Caption         =   "&View"
      Enabled         =   0   'False
   End
   Begin VB.Menu mnuTest 
      Caption         =   "&Test"
   End
End
Attribute VB_Name = "frmMain"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False


Option Explicit

Private Sub Form_Load()
    
    Dim lpRect As RECT
    Dim lRet As Long
        
    gsubScreenDimensions Me, 0, 0

    lpRect.Top = SCREEN_TOP
    lpRect.Left = SCREEN_LEFT + 800
    lpRect.Right = SCREEN_WIDTH + 800
    lpRect.Bottom = SCREEN_HEIGHT
    lRet = ClipCursor(lpRect)

End Sub


Private Sub Form_Unload(Cancel As Integer)
    
    Dim lpRect As RECT
    Dim lRet As Long
        
    lpRect.Top = 0
    lpRect.Left = 0
    lpRect.Right = (Screen.Width / Screen.TwipsPerPixelX) * 2
    lpRect.Bottom = Screen.Height / Screen.TwipsPerPixelY
    lRet = ClipCursor(lpRect)

End Sub


Private Sub mnuClose_Click()
    
    gsProfile = ""
    mnuDefaults.Enabled = False
    mnuSetup.Enabled = False
    mnuRun.Enabled = False
    mnuResults.Enabled = False
    
    If gsProfile = "" And gsResultsFile = "" Then
        mnuView.Enabled = False
    End If
    
    mnuNew.Enabled = True
    mnuOpen.Enabled = True
    mnuClose.Enabled = False

End Sub

Private Sub mnuCloseResults_Click()
    
    gsResultsFile = ""
    mnuDefaults.Enabled = False
    mnuSetup.Enabled = False
    mnuRun.Enabled = False
    mnuResults.Enabled = False
    
    If gsProfile = "" And gsResultsFile = "" Then
        mnuView.Enabled = False
    End If
    
    mnuNewResults.Enabled = True
    mnuOpenResults.Enabled = True
    mnuCloseResults.Enabled = False

End Sub


Private Sub mnuDefaults_Click()

    frmDefaults.Show 1

End Sub

Private Sub mnuExit_Click()

    
    Unload Me
    End
    
End Sub


Private Sub mnuNew_Click()
    
    On Error GoTo Error_mnuNew_Click
    
    diaFile.Filter = "Session DB (*.ts)|*.ts|All Files (*.*)|*.*"
    diaFile.DialogTitle = "New Change Detection Session Database"
    diaFile.DefaultExt = "ts"
    diaFile.InitDir = App.Path
    diaFile.CancelError = True
    diaFile.Flags = cdlOFNCreatePrompt + cdlOFNHideReadOnly + cdlOFNOverwritePrompt
    diaFile.ShowSave
    
    gsProfile = diaFile.FileName
    
    FileCopy gsProfile, CurDir & "\bak" & diaFile.FileTitle
    FileCopy CurDir & "\Session.tmp", gsProfile
    
    If gsProfile <> "" And gsResultsFile <> "" Then
        mnuDefaults.Enabled = True
    End If
    
    mnuView.Enabled = True
    
    mnuNew.Enabled = False
    mnuOpen.Enabled = False
    mnuClose.Enabled = True
    
    Exit Sub
    

Error_mnuNew_Click:

    If Err = 53 Then   'file error
        Resume Next
    Else
        Exit Sub
    End If

End Sub

Private Sub mnuNewResults_Click()
    
    On Error GoTo Error_mnuNewResults_Click
    
    diaResults.Filter = "Results DB (*.tr)|*.tr|All Files (*.*)|*.*"
    diaResults.DialogTitle = "New Change Detection Results Database"
    diaResults.DefaultExt = "tr"
    diaResults.InitDir = App.Path
    diaResults.CancelError = True
    diaResults.Flags = cdlOFNCreatePrompt + cdlOFNHideReadOnly + cdlOFNOverwritePrompt
    diaResults.ShowSave
    
    gsResultsFile = diaResults.FileName
    FileCopy gsResultsFile, CurDir & "\bak" & diaResults.FileTitle
    FileCopy CurDir & "\Results.tmp", gsResultsFile
    
    If gsProfile <> "" And gsResultsFile <> "" Then
        mnuDefaults.Enabled = True
    End If
    
    mnuView.Enabled = True
    
    mnuNewResults.Enabled = False
    mnuOpenResults.Enabled = False
    mnuCloseResults.Enabled = True
    
    Exit Sub
    

Error_mnuNewResults_Click:

    If Err = 53 Then  'file error
        Resume Next
    Else
        Exit Sub
    End If

End Sub


Private Sub mnuOpen_Click()
    
    On Error GoTo Error_mnuOpen_Click
    
    diaFile.Filter = "Session DB (*.ts)|*.ts|All Files (*.*)|*.*"
    diaFile.DialogTitle = "Open Change Detection Session Database"
    diaFile.DefaultExt = "ts"
    diaFile.InitDir = App.Path
    diaFile.CancelError = True
    diaFile.Flags = cdlOFNFileMustExist + cdlOFNHideReadOnly
    diaFile.ShowOpen
    
    gsProfile = diaFile.FileName
    FileCopy gsProfile, CurDir & "\bak" & diaFile.FileTitle
    
    If gsProfile <> "" And gsResultsFile <> "" Then
        mnuDefaults.Enabled = True
    End If
    
    mnuView.Enabled = True
    
    mnuNew.Enabled = False
    mnuOpen.Enabled = False
    mnuClose.Enabled = True
    
    Exit Sub
    
Error_mnuOpen_Click:

    If Err = 76 Then  'file error
        Resume Next
    Else
        Exit Sub
    End If
    
End Sub

Private Sub mnuOpenResults_Click()
    
    On Error GoTo Error_mnuOpenResults_Click
    
    diaResults.Filter = "Session DB (*.tr)|*.tr|All Files (*.*)|*.*"
    diaResults.DialogTitle = "Open Change Detection Session Database"
    diaResults.DefaultExt = "tr"
    diaResults.InitDir = App.Path
    diaResults.CancelError = True
    diaResults.Flags = cdlOFNFileMustExist + cdlOFNHideReadOnly
    diaResults.ShowOpen
    
    gsResultsFile = diaResults.FileName
    
    FileCopy gsResultsFile, CurDir & "\bak" & diaResults.FileTitle
    
    If gsProfile <> "" And gsResultsFile <> "" Then
        mnuDefaults.Enabled = True
    End If
    
    mnuView.Enabled = True
    
    mnuNewResults.Enabled = False
    mnuOpenResults.Enabled = False
    mnuCloseResults.Enabled = True
    
    Exit Sub
    
Error_mnuOpenResults_Click:

    If Err = 53 Then  'file error
        Resume Next
    Else
        Exit Sub
    End If

End Sub


Private Sub mnuResults_Click()

    frmResults.Show 1
    
End Sub

Private Sub mnuRun_Click()

    frmControl.Show
    
End Sub

Private Sub mnuSetup_Click()

    frmSetup.Show 1
    
End Sub


Private Sub mnuTest_Click()

    frmTest.Show 1
    
End Sub


Private Sub mnuView_Click()

    frmViewData.Show 1
    
End Sub


