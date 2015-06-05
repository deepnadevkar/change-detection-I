VERSION 5.00
Object = "{8D650141-6025-11D1-BC40-0000C042AEC0}#3.0#0"; "ssdw3b32.ocx"
Begin VB.Form frmResults 
   BackColor       =   &H00404040&
   Caption         =   "Results"
   ClientHeight    =   6255
   ClientLeft      =   555
   ClientTop       =   1020
   ClientWidth     =   9165
   ControlBox      =   0   'False
   LinkTopic       =   "Form1"
   MaxButton       =   0   'False
   MinButton       =   0   'False
   PaletteMode     =   1  'UseZOrder
   ScaleHeight     =   6255
   ScaleWidth      =   9165
   Begin VB.Data datStimResponses 
      Caption         =   "StimResponses"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   6750
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   5160
      Visible         =   0   'False
      Width           =   2385
   End
   Begin VB.Data datDefaults 
      Caption         =   "Defaults"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   2760
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   "SessionInfo"
      Top             =   5100
      Visible         =   0   'False
      Width           =   2385
   End
   Begin VB.Data datSessionResults 
      Caption         =   "Session Results"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   2760
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   "SessionResults"
      Top             =   4440
      Visible         =   0   'False
      Width           =   2385
   End
   Begin VB.Data datTrialResults 
      Caption         =   "Trial Results"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   2760
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   "TrialResults"
      Top             =   4770
      Visible         =   0   'False
      Width           =   2385
   End
   Begin VB.CheckBox chkVariableListLength 
      BackColor       =   &H00404040&
      DataField       =   "VariableListLength"
      DataSource      =   "datSessionResults"
      Enabled         =   0   'False
      Height          =   255
      Left            =   2490
      TabIndex        =   25
      Top             =   1110
      Width           =   315
   End
   Begin VB.CheckBox chkAbortCorrectionEnabled 
      BackColor       =   &H00404040&
      DataField       =   "AbortCorrectionEnabled"
      DataSource      =   "datSessionResults"
      Enabled         =   0   'False
      Height          =   255
      Left            =   2490
      TabIndex        =   1
      Top             =   1440
      Width           =   315
   End
   Begin VB.CheckBox chkIncorrectCorrectionEnabled 
      BackColor       =   &H00404040&
      DataField       =   "IncorrectCorrectionEnabled"
      DataSource      =   "datSessionResults"
      Enabled         =   0   'False
      Height          =   255
      Left            =   2490
      TabIndex        =   0
      Top             =   2100
      Width           =   285
   End
   Begin SSDataWidgets_B.SSDBGrid gridStimResponses 
      Bindings        =   "RESULTS.frx":0000
      Height          =   3525
      Left            =   7200
      TabIndex        =   56
      Top             =   2940
      Width           =   2235
      _Version        =   196617
      BeginProperty HeadFont {0BE35203-8F91-11CE-9DE3-00AA004BB851} 
         Name            =   "MS Sans Serif"
         Size            =   8.25
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      AllowUpdate     =   0   'False
      MultiLine       =   0   'False
      AllowRowSizing  =   0   'False
      AllowGroupSizing=   0   'False
      AllowColumnSizing=   0   'False
      AllowGroupMoving=   0   'False
      AllowColumnMoving=   0
      AllowGroupSwapping=   0   'False
      AllowColumnSwapping=   0
      AllowGroupShrinking=   0   'False
      AllowColumnShrinking=   0   'False
      AllowDragDrop   =   0   'False
      SelectTypeRow   =   1
      BalloonHelp     =   0   'False
      ForeColorEven   =   0
      BackColorOdd    =   16777215
      RowHeight       =   423
      Columns.Count   =   2
      Columns(0).Width=   847
      Columns(0).Caption=   "Pos"
      Columns(0).Alignment=   1
      Columns(0).CaptionAlignment=   1
      Columns(0).DataField=   "ListPosition"
      Columns(0).DataType=   3
      Columns(0).FieldLen=   256
      Columns(1).Width=   2117
      Columns(1).Caption=   "RespTime"
      Columns(1).Alignment=   1
      Columns(1).CaptionAlignment=   1
      Columns(1).DataField=   "ResponseTime"
      Columns(1).DataType=   3
      Columns(1).FieldLen=   256
      UseDefaults     =   0   'False
      _ExtentX        =   3942
      _ExtentY        =   6218
      _StockProps     =   79
      Caption         =   "Stim Responses"
      BeginProperty Font {0BE35203-8F91-11CE-9DE3-00AA004BB851} 
         Name            =   "MS Sans Serif"
         Size            =   8.25
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      BeginProperty PageFooterFont {0BE35203-8F91-11CE-9DE3-00AA004BB851} 
         Name            =   "MS Sans Serif"
         Size            =   8.25
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      BeginProperty PageHeaderFont {0BE35203-8F91-11CE-9DE3-00AA004BB851} 
         Name            =   "MS Sans Serif"
         Size            =   8.25
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
   End
   Begin SSDataWidgets_B.SSDBGrid gridTrialResults 
      Bindings        =   "RESULTS.frx":001F
      Height          =   3555
      Left            =   300
      TabIndex        =   55
      Top             =   2940
      Width           =   6855
      _Version        =   196617
      BeginProperty HeadFont {0BE35203-8F91-11CE-9DE3-00AA004BB851} 
         Name            =   "MS Sans Serif"
         Size            =   8.25
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      GroupHeadLines  =   2
      AllowUpdate     =   0   'False
      AllowRowSizing  =   0   'False
      AllowGroupSizing=   0   'False
      AllowColumnSizing=   0   'False
      AllowGroupMoving=   0   'False
      AllowColumnMoving=   0
      AllowGroupSwapping=   0   'False
      AllowColumnSwapping=   0
      AllowGroupShrinking=   0   'False
      AllowColumnShrinking=   0   'False
      AllowDragDrop   =   0   'False
      SelectTypeRow   =   1
      ForeColorEven   =   0
      BackColorOdd    =   16777215
      RowHeight       =   423
      Columns.Count   =   7
      Columns(0).Width=   926
      Columns(0).Caption=   "Trial"
      Columns(0).Alignment=   1
      Columns(0).CaptionAlignment=   2
      Columns(0).DataField=   "TrialNum"
      Columns(0).DataType=   3
      Columns(0).FieldLen=   256
      Columns(1).Width=   1138
      Columns(1).Caption=   "Type"
      Columns(1).Alignment=   1
      Columns(1).CaptionAlignment=   2
      Columns(1).DataField=   "TrialType"
      Columns(1).DataType=   8
      Columns(1).FieldLen=   256
      Columns(2).Width=   1720
      Columns(2).Caption=   "Response"
      Columns(2).Alignment=   1
      Columns(2).CaptionAlignment=   2
      Columns(2).DataField=   "Response"
      Columns(2).DataType=   8
      Columns(2).FieldLen=   256
      Columns(3).Width=   1773
      Columns(3).Caption=   "Resp Time"
      Columns(3).Alignment=   1
      Columns(3).CaptionAlignment=   2
      Columns(3).DataField=   "ResponseTime"
      Columns(3).DataType=   3
      Columns(3).FieldLen=   256
      Columns(4).Width=   1826
      Columns(4).Caption=   "Cum % Corr"
      Columns(4).Alignment=   1
      Columns(4).CaptionAlignment=   2
      Columns(4).DataField=   "CumPercentCorrect"
      Columns(4).DataType=   4
      Columns(4).FieldLen=   256
      Columns(5).Width=   1799
      Columns(5).Caption=   "AbortCorr"
      Columns(5).Alignment=   1
      Columns(5).CaptionAlignment=   2
      Columns(5).DataField=   "AbortCorrections"
      Columns(5).DataType=   3
      Columns(5).FieldLen=   256
      Columns(6).Width=   1879
      Columns(6).Caption=   "IncorrCorr"
      Columns(6).Alignment=   1
      Columns(6).CaptionAlignment=   2
      Columns(6).DataField=   "IncorrectCorrections"
      Columns(6).DataType=   3
      Columns(6).FieldLen=   256
      UseDefaults     =   0   'False
      _ExtentX        =   12091
      _ExtentY        =   6271
      _StockProps     =   79
      Caption         =   "Trial Results"
      BeginProperty Font {0BE35203-8F91-11CE-9DE3-00AA004BB851} 
         Name            =   "MS Sans Serif"
         Size            =   8.25
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      BeginProperty PageFooterFont {0BE35203-8F91-11CE-9DE3-00AA004BB851} 
         Name            =   "MS Sans Serif"
         Size            =   8.25
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      BeginProperty PageHeaderFont {0BE35203-8F91-11CE-9DE3-00AA004BB851} 
         Name            =   "MS Sans Serif"
         Size            =   8.25
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
   End
   Begin VB.Label Label39 
      BackColor       =   &H00404040&
      Caption         =   "Diff Incorrect:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   6750
      TabIndex        =   54
      Top             =   2460
      Width           =   1305
   End
   Begin VB.Label Label38 
      BackColor       =   &H00404040&
      Caption         =   "Same Incorrect:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   6750
      TabIndex        =   53
      Top             =   2130
      Width           =   1305
   End
   Begin VB.Label lblSameIncorrect 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "SameIncorrect"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   8370
      TabIndex        =   52
      Top             =   2100
      Width           =   885
   End
   Begin VB.Label lblDiffIncorrect 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "DiffIncorrect"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   8370
      TabIndex        =   51
      Top             =   2430
      Width           =   885
   End
   Begin VB.Label Label35 
      BackColor       =   &H00404040&
      Caption         =   "Diff Correct:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   6750
      TabIndex        =   50
      Top             =   1800
      Width           =   1305
   End
   Begin VB.Label Label34 
      BackColor       =   &H00404040&
      Caption         =   "Same Correct:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   6750
      TabIndex        =   49
      Top             =   1470
      Width           =   1305
   End
   Begin VB.Label Label33 
      BackColor       =   &H00404040&
      Caption         =   "Percent Correct:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   6750
      TabIndex        =   48
      Top             =   1140
      Width           =   1305
   End
   Begin VB.Label Label32 
      BackColor       =   &H00404040&
      Caption         =   "Incorrect Corrections:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   6750
      TabIndex        =   47
      Top             =   810
      Width           =   1575
   End
   Begin VB.Label Label20 
      BackColor       =   &H00404040&
      Caption         =   "Abort Corrections:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   6750
      TabIndex        =   46
      Top             =   480
      Width           =   1305
   End
   Begin VB.Label lblAbortCorrections 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "TotalAbortCorrections"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   8370
      TabIndex        =   45
      Top             =   450
      Width           =   885
   End
   Begin VB.Label lblIncorrectCorrections 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "TotalIncorrectCorrections"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   8370
      TabIndex        =   44
      Top             =   780
      Width           =   885
   End
   Begin VB.Label lblPercentCorrect 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "PercentCorrect"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   8370
      TabIndex        =   43
      Top             =   1110
      Width           =   885
   End
   Begin VB.Label lblSameCorrect 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "SameCorrect"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   8370
      TabIndex        =   42
      Top             =   1440
      Width           =   885
   End
   Begin VB.Label lblDiffCorrect 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "DiffCorrect"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   8370
      TabIndex        =   41
      Top             =   1770
      Width           =   885
   End
   Begin VB.Label Label11 
      BackColor       =   &H00404040&
      Caption         =   "Aborts:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   6750
      TabIndex        =   40
      Top             =   150
      Width           =   1065
   End
   Begin VB.Label lblAborts 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "TotalAborts"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   8370
      TabIndex        =   39
      Top             =   120
      Width           =   885
   End
   Begin VB.Label Label12 
      BackColor       =   &H00404040&
      Caption         =   "End:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   3570
      TabIndex        =   38
      Top             =   810
      Width           =   525
   End
   Begin VB.Label lblEndTime 
      Alignment       =   2  'Center
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "EndTime"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   4200
      TabIndex        =   37
      Top             =   780
      Width           =   1725
   End
   Begin VB.Label Label10 
      BackColor       =   &H00404040&
      Caption         =   "Start:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   3570
      TabIndex        =   36
      Top             =   480
      Width           =   525
   End
   Begin VB.Label lblStartTime 
      Alignment       =   2  'Center
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "StartTime"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   4200
      TabIndex        =   35
      Top             =   450
      Width           =   1725
   End
   Begin VB.Label Label2 
      BackColor       =   &H00404040&
      Caption         =   "Subject:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   330
      TabIndex        =   34
      Top             =   480
      Width           =   675
   End
   Begin VB.Label Label3 
      BackColor       =   &H00404040&
      Caption         =   "Session ID:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   330
      TabIndex        =   33
      Top             =   150
      Width           =   915
   End
   Begin VB.Label Label4 
      BackColor       =   &H00404040&
      Caption         =   "Number of Trials:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   330
      TabIndex        =   32
      Top             =   810
      Width           =   1275
   End
   Begin VB.Label Label1 
      BackColor       =   &H00404040&
      Caption         =   "Date:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   3570
      TabIndex        =   31
      Top             =   150
      Width           =   525
   End
   Begin VB.Label lblNumOfTrials 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "NumOfTrials"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   1800
      TabIndex        =   30
      Top             =   780
      Width           =   885
   End
   Begin VB.Label lblSubject 
      Alignment       =   2  'Center
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "Subject"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   1800
      TabIndex        =   29
      Top             =   450
      Width           =   885
   End
   Begin VB.Label lblSessionID 
      Alignment       =   2  'Center
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "SessionID"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   1800
      TabIndex        =   28
      Top             =   120
      Width           =   885
   End
   Begin VB.Label lblDate 
      Alignment       =   2  'Center
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "Date"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   4200
      TabIndex        =   27
      Top             =   120
      Width           =   1725
   End
   Begin VB.Label Label31 
      BackColor       =   &H00404040&
      Caption         =   "Variable ListLength:"
      ForeColor       =   &H00C0C0C0&
      Height          =   225
      Left            =   330
      TabIndex        =   26
      Top             =   1140
      Width           =   1515
   End
   Begin VB.Label Label30 
      BackColor       =   &H00404040&
      Caption         =   "%"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   6030
      TabIndex        =   24
      Top             =   2430
      Width           =   345
   End
   Begin VB.Label Label29 
      BackColor       =   &H00404040&
      Caption         =   "%"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   6030
      TabIndex        =   23
      Top             =   2100
      Width           =   345
   End
   Begin VB.Label Label28 
      BackColor       =   &H00404040&
      Caption         =   "%"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   6030
      TabIndex        =   22
      Top             =   1770
      Width           =   345
   End
   Begin VB.Label Label27 
      BackColor       =   &H00404040&
      Caption         =   "%"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   6030
      TabIndex        =   21
      Top             =   1440
      Width           =   345
   End
   Begin VB.Label Label26 
      BackColor       =   &H00404040&
      Caption         =   "%"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   6030
      TabIndex        =   20
      Top             =   1110
      Width           =   345
   End
   Begin VB.Label Label25 
      BackColor       =   &H00404040&
      Caption         =   "Diff Probability:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   3570
      TabIndex        =   19
      Top             =   2460
      Width           =   1305
   End
   Begin VB.Label Label24 
      BackColor       =   &H00404040&
      Caption         =   "Same Probability:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   3570
      TabIndex        =   18
      Top             =   2130
      Width           =   1305
   End
   Begin VB.Label Label23 
      BackColor       =   &H00404040&
      Caption         =   "Reinf Probability:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   3570
      TabIndex        =   17
      Top             =   1800
      Width           =   1305
   End
   Begin VB.Label Label22 
      BackColor       =   &H00404040&
      Caption         =   "Pellet Probability:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   3570
      TabIndex        =   16
      Top             =   1470
      Width           =   1305
   End
   Begin VB.Label Label21 
      BackColor       =   &H00404040&
      Caption         =   "Juice Probability:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   3570
      TabIndex        =   15
      Top             =   1140
      Width           =   1305
   End
   Begin VB.Label lblJuiceProbability 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "JuiceProbability"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   5040
      TabIndex        =   14
      Top             =   1110
      Width           =   885
   End
   Begin VB.Label lblPelletProbability 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "PelletProbability"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   5040
      TabIndex        =   13
      Top             =   1440
      Width           =   885
   End
   Begin VB.Label lblReinfProbability 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "ReinfProbability"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   5040
      TabIndex        =   12
      Top             =   1770
      Width           =   885
   End
   Begin VB.Label lblSameProbability 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "SameProbability"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   5040
      TabIndex        =   11
      Top             =   2100
      Width           =   885
   End
   Begin VB.Label lblDiffProbability 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "DiffProbability"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   5040
      TabIndex        =   10
      Top             =   2430
      Width           =   885
   End
   Begin VB.Label Label5 
      BackColor       =   &H00404040&
      Caption         =   "Abort Corrections Enabled:"
      ForeColor       =   &H00C0C0C0&
      Height          =   225
      Left            =   330
      TabIndex        =   9
      Top             =   1470
      Width           =   2025
   End
   Begin VB.Label Label6 
      BackColor       =   &H00404040&
      Caption         =   "Incorrect Corrections Enabled:"
      ForeColor       =   &H00C0C0C0&
      Height          =   225
      Left            =   330
      TabIndex        =   8
      Top             =   2130
      Width           =   2205
   End
   Begin VB.Label Label7 
      BackColor       =   &H00404040&
      Caption         =   "Abort LTO:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   330
      TabIndex        =   7
      Top             =   1800
      Width           =   1065
   End
   Begin VB.Label Label8 
      BackColor       =   &H00404040&
      Caption         =   "Incorrect LTO:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   330
      TabIndex        =   6
      Top             =   2460
      Width           =   1065
   End
   Begin VB.Label Label15 
      BackColor       =   &H00404040&
      Caption         =   "sec."
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   2790
      TabIndex        =   5
      Top             =   1800
      Width           =   405
   End
   Begin VB.Label Label16 
      BackColor       =   &H00404040&
      Caption         =   "sec."
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   2790
      TabIndex        =   4
      Top             =   2460
      Width           =   405
   End
   Begin VB.Label lblAbortLTO 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "DefaultAbortLTO"
      DataSource      =   "datDefaults"
      Height          =   285
      Left            =   1800
      TabIndex        =   3
      Top             =   1770
      Width           =   885
   End
   Begin VB.Label lblIncorrectLTO 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "DefaultIncorrectLTO"
      DataSource      =   "datDefaults"
      Height          =   285
      Left            =   1800
      TabIndex        =   2
      Top             =   2430
      Width           =   885
   End
   Begin VB.Menu mnuFile 
      Caption         =   "&File"
      Begin VB.Menu mnuPrint 
         Caption         =   "&Print"
      End
      Begin VB.Menu mnuSep1 
         Caption         =   "-"
      End
      Begin VB.Menu mnuClose 
         Caption         =   "&Close"
      End
   End
End
Attribute VB_Name = "frmResults"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False
Option Explicit

Private Sub subPrintResults()

    Dim bCheck1 As Boolean
    Dim bCheck2 As Boolean
    Dim bCheck3 As Boolean
    Dim sField1 As String
    Dim sField2 As String
    Dim sField3 As String
    Dim sField4 As String
    Dim sField5 As String
    Dim sField6 As String
    Dim sField7 As String
    
    sField1 = Space(6)
    sField2 = sField1
    sField3 = sField1
    sField4 = sField1
    sField5 = sField1
    sField6 = sField1
    sField7 = sField1
    
    If chkVariableListLength = 0 Then
        bCheck1 = False
    Else
        bCheck1 = True
    End If
    
    If chkAbortCorrectionEnabled = 0 Then
        bCheck2 = False
    Else
        bCheck2 = True
    End If
    
    If chkIncorrectCorrectionEnabled = 0 Then
        bCheck3 = False
    Else
        bCheck3 = True
    End If
    
    If datSessionResults.Recordset.RecordCount <> 0 Then

        Printer.Font.name = "Courier"
        
        Printer.Print
        Printer.Print
        Printer.Print
        Printer.Print Tab(10); "SPR Training"
        Printer.Print Tab(10); "SPR Results: " & gsResultsFile
        Printer.Print Tab(10); "Phase " & datDefaults.Recordset("Phase")
        
        Printer.Print
        Printer.Print
        Printer.Print
        Printer.Print Tab(10); "SessionID:"; Tab(35); lblSessionID
        Printer.Print Tab(10); "Subject:"; Tab(35); lblSubject
        Printer.Print Tab(10); "NumOfTrials:"; Tab(35); lblNumOfTrials
        Printer.Print Tab(10); "Variable LengthList:"; Tab(35); bCheck1
        Printer.Print Tab(10); "Abort Corr Enabled:"; Tab(35); bCheck2
        Printer.Print Tab(10); "Abort LTO:"; Tab(35); lblAbortLTO
        Printer.Print Tab(10); "Incorrect Corr Enabled:"; Tab(35); bCheck3
        Printer.Print Tab(10); "Incorrect LTO:"; Tab(35); lblIncorrectLTO
        Printer.Print Tab(10); "Date:"; Tab(35); lblDate
        Printer.Print Tab(10); "Start:"; Tab(35); lblStartTime
        Printer.Print Tab(10); "End:"; Tab(35); lblEndTime
        Printer.Print Tab(10); "Juice Probability:"; Tab(35); lblJuiceProbability
        Printer.Print Tab(10); "Pellet Probability:"; Tab(35); lblPelletProbability
        Printer.Print Tab(10); "Reinf Probability:"; Tab(35); lblReinfProbability
        Printer.Print Tab(10); "Same Probability:"; Tab(35); lblSameProbability
        Printer.Print Tab(10); "Diff Probability:"; Tab(35); lblDiffProbability
        Printer.Print Tab(10); "Aborts:"; Tab(35); lblAborts
        Printer.Print Tab(10); "Abort Corr:"; Tab(35); lblAbortCorrections
        Printer.Print Tab(10); "Incorrect Corr:"; Tab(35); lblIncorrectCorrections
        Printer.Print Tab(10); "Percent Correct:"; Tab(35); lblPercentCorrect
        Printer.Print Tab(10); "Same Correct:"; Tab(35); lblSameCorrect
        Printer.Print Tab(10); "Diff Correct:"; Tab(35); lblDiffCorrect
        Printer.Print Tab(10); "Same Incorrect:"; Tab(35); lblSameIncorrect
        Printer.Print Tab(10); "Diff Incorrect:"; Tab(35); lblDiffIncorrect
    
    End If
    
    Printer.Print
    Printer.Print
    Printer.Print
    Printer.Print
    
    'heading
    Printer.Print Tab(10); "       "; Tab(20); "      "; Tab(30); "      "; Tab(40); "  Resp"; Tab(50); "  Percent"; Tab(60); "  Abort"; Tab(70); "  Incorr"
    Printer.Print Tab(10); "  Trial"; Tab(20); "  Type"; Tab(30); "  Resp"; Tab(40); "  Time"; Tab(50); "  Correct"; Tab(60); "  Corr "; Tab(70); "  Corr  "
    Printer.Print
    
    If datTrialResults.Recordset.RecordCount <> 0 Then
        datTrialResults.Recordset.MoveFirst
    
        While Not datTrialResults.Recordset.EOF
            RSet sField1 = datTrialResults.Recordset("TrialNum")
            RSet sField2 = datTrialResults.Recordset("TrialType")
            RSet sField3 = datTrialResults.Recordset("Response")
            RSet sField4 = Format(datTrialResults.Recordset("ResponseTime") / 1000, "##0.0")
            RSet sField5 = datTrialResults.Recordset("CumPercentCorrect")
            RSet sField6 = datTrialResults.Recordset("AbortCorrections")
            RSet sField7 = datTrialResults.Recordset("IncorrectCorrections")
            Printer.Print Tab(10); sField1; _
                        Tab(20); sField2; _
                        Tab(30); sField3; _
                        Tab(40); sField4; _
                        Tab(50); sField5; _
                        Tab(60); sField6; _
                        Tab(70); sField7
        
            datTrialResults.Recordset.MoveNext
            
        Wend
    End If
    
    Printer.EndDoc
    
End Sub

Private Sub Form_Load()
    
    gsubScreenDimensions Me, 0, 0
    
    datSessionResults.DatabaseName = gsResultsFile
    datSessionResults.RecordSource = "SessionResults"
    datSessionResults.Refresh
    
    datTrialResults.DatabaseName = gsResultsFile
    datTrialResults.RecordSource = "TrialResults"
    datTrialResults.Refresh
    
    datDefaults.DatabaseName = gsProfile
    datDefaults.RecordSource = "SessionInfo"
    datDefaults.Refresh
    gridTrialResults.Caption = "Trial Results: Phase " & datDefaults.Recordset("Phase")
    
    If datSessionResults.Recordset.RecordCount = 0 And _
        datTrialResults.Recordset.RecordCount = 0 And _
        datDefaults.Recordset.RecordCount = 0 Then
       
        mnuPrint.Enabled = False
       
    Else
    
        mnuPrint.Enabled = True
        
    End If

End Sub


Private Sub gridTrialResults_Click()

    datStimResponses.DatabaseName = gsResultsFile
    datStimResponses.RecordSource = "SELECT * FROM StimResponses WHERE TrialNum = " & datTrialResults.Recordset("TrialNum")
    datStimResponses.Refresh
    
End Sub

Private Sub mnuClose_Click()

    Unload Me
    
End Sub


Private Sub mnuPrint_Click()

    subPrintResults
    
End Sub


