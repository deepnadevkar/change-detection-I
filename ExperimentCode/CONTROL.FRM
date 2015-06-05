VERSION 5.00
Object = "{02B97C13-ED1E-11CD-A08B-00AA00575482}#1.0#0"; "mhiinp32.ocx"
Object = "{3AE5AE83-A6DA-101B-9313-00AA00575482}#1.0#0"; "Mhfram32.ocx"
Object = "{C30E0AD3-B1BA-11CE-ABC6-F5B2E79D9E3F}#6.0#0"; "lead32.ocx"
Object = "{B08ECCE3-EE86-11CE-9577-0020AF039CA3}#1.0#0"; "mhwave32.ocx"
Object = "{4DE9E2A3-150F-11CF-8FBF-444553540000}#4.0#0"; "DlxOCX32.ocx"
Begin VB.Form frmControl 
   BackColor       =   &H00404040&
   BorderStyle     =   1  'Fixed Single
   Caption         =   "Control Panel"
   ClientHeight    =   7620
   ClientLeft      =   810
   ClientTop       =   1215
   ClientWidth     =   10245
   ControlBox      =   0   'False
   LinkTopic       =   "Form1"
   MaxButton       =   0   'False
   MinButton       =   0   'False
   PaletteMode     =   1  'UseZOrder
   ScaleHeight     =   7620
   ScaleWidth      =   10245
   Begin MhinintLibCtl.MhIntInput intJuiceDelay 
      DataField       =   "JuiceTime"
      DataSource      =   "datDefaults"
      Height          =   315
      Left            =   8280
      TabIndex        =   85
      Top             =   2220
      Visible         =   0   'False
      Width           =   585
      _Version        =   65536
      _ExtentX        =   1032
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Text            =   "1000"
      Max             =   30000
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intPelletDelay 
      DataField       =   "PelletTime"
      DataSource      =   "datDefaults"
      Height          =   315
      Left            =   8280
      TabIndex        =   83
      Top             =   2550
      Width           =   585
      _Version        =   65536
      _ExtentX        =   1032
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Max             =   30000
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhwaveLib.MhWave waveBeep 
      Height          =   420
      Left            =   2790
      TabIndex        =   81
      Top             =   1050
      Visible         =   0   'False
      Width           =   420
      _Version        =   65536
      _ExtentX        =   741
      _ExtentY        =   741
      _StockProps     =   0
      FilenameRead    =   "CONTROL.frx":0000
   End
   Begin VB.Timer tmrAfterReady 
      Enabled         =   0   'False
      Left            =   840
      Top             =   6720
   End
   Begin DlsrLib.DriverLINXSR DriverLINXSR1 
      Left            =   2760
      Top             =   1680
      _Version        =   262144
      _ExtentX        =   741
      _ExtentY        =   741
      _StockProps     =   64
   End
   Begin VB.Timer tmrHouseLight1 
      Left            =   2880
      Top             =   2880
   End
   Begin VB.Timer tmrHouseLight 
      Left            =   2880
      Top             =   2280
   End
   Begin VB.Timer tmrCorrectItemDuration 
      Enabled         =   0   'False
      Left            =   6600
      Tag             =   "CorrectItemDuration"
      Top             =   5160
   End
   Begin VB.Timer tmrSampleDuration 
      Enabled         =   0   'False
      Left            =   6720
      Tag             =   "SampleDuration"
      Top             =   4680
   End
   Begin VB.CheckBox Check3 
      BackColor       =   &H00404040&
      DataField       =   "ProbeAbortEnabled"
      DataSource      =   "datSessionInfo"
      Enabled         =   0   'False
      Height          =   255
      Left            =   8700
      TabIndex        =   109
      Top             =   7050
      Visible         =   0   'False
      Width           =   285
   End
   Begin VB.Timer tmrProbeAbort 
      Enabled         =   0   'False
      Left            =   6810
      Tag             =   "ProbeAbort"
      Top             =   4140
   End
   Begin VB.CheckBox Check2 
      BackColor       =   &H00404040&
      DataField       =   "TouchAbortEnabled"
      DataSource      =   "datSessionInfo"
      Enabled         =   0   'False
      Height          =   255
      Left            =   5790
      TabIndex        =   107
      Top             =   6780
      Visible         =   0   'False
      Width           =   285
   End
   Begin VB.CheckBox Check1 
      BackColor       =   &H00404040&
      DataField       =   "NoTouchProbes"
      DataSource      =   "datSessionInfo"
      Enabled         =   0   'False
      Height          =   255
      Left            =   5940
      TabIndex        =   105
      Top             =   6690
      Visible         =   0   'False
      Width           =   255
   End
   Begin VB.Data datStimResponses 
      Caption         =   "StimResponses"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   345
      Left            =   3990
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   2820
      Visible         =   0   'False
      Width           =   2475
   End
   Begin VB.Timer tmrApparatusDelay 
      Enabled         =   0   'False
      Left            =   8880
      Top             =   1050
   End
   Begin VB.Data datRecordCounter 
      Caption         =   "RecordCounter"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   345
      Left            =   4260
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   630
      Visible         =   0   'False
      Width           =   2655
   End
   Begin VB.Timer tmrMetraByte 
      Enabled         =   0   'False
      Left            =   6780
      Top             =   870
   End
   Begin VB.Timer tmrInit 
      Enabled         =   0   'False
      Left            =   6780
      Top             =   420
   End
   Begin VB.Timer tmrReadTouch 
      Enabled         =   0   'False
      Interval        =   100
      Left            =   7110
      Top             =   0
   End
   Begin VB.Timer tmrClock 
      Interval        =   1000
      Left            =   2760
      Top             =   480
   End
   Begin VB.Timer tmrITI 
      Enabled         =   0   'False
      Left            =   6780
      Tag             =   "ITI"
      Top             =   3690
   End
   Begin VB.Timer tmrIncorrectLTO 
      Enabled         =   0   'False
      Left            =   6840
      Tag             =   "IncorrectLTO"
      Top             =   3270
   End
   Begin VB.Timer tmrAbortLTO 
      Enabled         =   0   'False
      Left            =   6840
      Tag             =   "AbortLTO"
      Top             =   2850
   End
   Begin VB.Timer tmrProbeView 
      Enabled         =   0   'False
      Left            =   6420
      Tag             =   "ProbeView"
      Top             =   2610
   End
   Begin VB.Timer tmrProbeDelay 
      Enabled         =   0   'False
      Left            =   6450
      Tag             =   "ProbeDelay"
      Top             =   2160
   End
   Begin VB.Timer tmrISI 
      Enabled         =   0   'False
      Left            =   6780
      Tag             =   "ISI"
      Top             =   1680
   End
   Begin VB.Timer tmrStimuliView 
      Enabled         =   0   'False
      Left            =   6690
      Tag             =   "StimuliView"
      Top             =   1290
   End
   Begin VB.Data datTrialResults 
      Caption         =   "TrialResults"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   345
      Left            =   4710
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   930
      Visible         =   0   'False
      Width           =   2085
   End
   Begin VB.Data datSessionInfo 
      Caption         =   "SessionInfo"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   345
      Left            =   4560
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   1560
      Visible         =   0   'False
      Width           =   2085
   End
   Begin VB.Data datSessionResults 
      Caption         =   "SessionResults"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   345
      Left            =   4710
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   1260
      Visible         =   0   'False
      Width           =   2085
   End
   Begin VB.Data datProbeInfo 
      Caption         =   "Probe"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   345
      Left            =   4710
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   2460
      Visible         =   0   'False
      Width           =   2085
   End
   Begin VB.Data datStimuliInfo 
      Caption         =   "Stimuli"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   345
      Left            =   4410
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   2250
      Visible         =   0   'False
      Width           =   2085
   End
   Begin VB.Data datTrialInfo 
      Caption         =   "TrialInfo"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   345
      Left            =   4440
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   1890
      Visible         =   0   'False
      Width           =   2085
   End
   Begin VB.CheckBox chkIncorrectCorrectionEnabled 
      BackColor       =   &H00404040&
      DataField       =   "IncorrectCorrectionsEnabled"
      DataSource      =   "datSessionInfo"
      Enabled         =   0   'False
      Height          =   255
      Left            =   2340
      TabIndex        =   2
      Top             =   5700
      Width           =   285
   End
   Begin VB.CheckBox chkAbortCorrectionEnabled 
      BackColor       =   &H00404040&
      DataField       =   "AbortCorrectionsEnabled"
      DataSource      =   "datSessionInfo"
      Enabled         =   0   'False
      Height          =   255
      Left            =   5820
      TabIndex        =   1
      Top             =   6900
      Visible         =   0   'False
      Width           =   315
   End
   Begin VB.CheckBox chkVariableProbeDelay 
      BackColor       =   &H00404040&
      DataField       =   "VariableProbeDelay"
      DataSource      =   "datSessionInfo"
      Enabled         =   0   'False
      Height          =   255
      Left            =   8700
      TabIndex        =   0
      Top             =   6840
      Visible         =   0   'False
      Width           =   315
   End
   Begin Mh3dfrmLibCtl.Mh3dFrame Mh3dFrame1 
      Height          =   2685
      Left            =   3360
      TabIndex        =   64
      Top             =   3840
      Width           =   3495
      _Version        =   65536
      _ExtentX        =   6165
      _ExtentY        =   4736
      _StockProps     =   79
      Caption         =   "Mh3dFrame1"
      BeginProperty Font {0BE35203-8F91-11CE-9DE3-00AA004BB851} 
         Name            =   "MS Sans Serif"
         Size            =   8.25
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Alignment       =   0
      AutoSize        =   0   'False
      BevelSize       =   4
      BevelStyle      =   0
      BorderColor     =   -2147483642
      BorderStyle     =   0
      FillColor       =   0
      FontStyle       =   0
      FontTransparent =   0   'False
      LightColor      =   -2147483643
      ShadowColor     =   -2147483632
      TextColor       =   -2147483640
      WallPaper       =   0
      NoPrefix        =   0   'False
      FormatString    =   ""
      Caption         =   "Mh3dFrame1"
      Picture         =   "CONTROL.frx":B820
      Begin VB.Label lblDiffTrialDiffResp 
         Alignment       =   2  'Center
         BackColor       =   &H000080FF&
         Caption         =   "0"
         DataField       =   "DiffCorrect"
         DataSource      =   "datSessionResults"
         BeginProperty Font 
            Name            =   "Arial"
            Size            =   18
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   435
         Left            =   2010
         TabIndex        =   95
         Top             =   1770
         Width           =   705
      End
      Begin VB.Label lblSameTrialDiffResp 
         Alignment       =   2  'Center
         BackColor       =   &H000080FF&
         Caption         =   "0"
         DataField       =   "SameIncorrect"
         DataSource      =   "datSessionResults"
         BeginProperty Font 
            Name            =   "Arial"
            Size            =   18
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   435
         Left            =   2010
         TabIndex        =   93
         Top             =   900
         Width           =   705
      End
      Begin VB.Label lblDiffTrialSameResp 
         Alignment       =   2  'Center
         BackColor       =   &H000080FF&
         Caption         =   "0"
         DataField       =   "DiffIncorrect"
         DataSource      =   "datSessionResults"
         BeginProperty Font 
            Name            =   "Arial"
            Size            =   18
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   435
         Left            =   1110
         TabIndex        =   94
         Top             =   1770
         Width           =   705
      End
      Begin VB.Label Label36 
         BackColor       =   &H000080FF&
         BorderStyle     =   1  'Fixed Single
         Height          =   855
         Left            =   1920
         TabIndex        =   98
         Top             =   1560
         Width           =   885
      End
      Begin VB.Label Label34 
         BackColor       =   &H000080FF&
         BorderStyle     =   1  'Fixed Single
         Height          =   855
         Left            =   1020
         TabIndex        =   97
         Top             =   1560
         Width           =   885
      End
      Begin VB.Label Label33 
         BackColor       =   &H000080FF&
         BorderStyle     =   1  'Fixed Single
         Height          =   855
         Left            =   1920
         TabIndex        =   96
         Top             =   690
         Width           =   885
      End
      Begin VB.Label lblSameTrialSameResp 
         Alignment       =   2  'Center
         BackColor       =   &H000080FF&
         Caption         =   "0"
         DataField       =   "SameCorrect"
         DataSource      =   "datSessionResults"
         BeginProperty Font 
            Name            =   "Arial"
            Size            =   18
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   435
         Left            =   1110
         TabIndex        =   92
         Top             =   900
         Width           =   705
      End
      Begin VB.Label Label32 
         BackColor       =   &H000080FF&
         BorderStyle     =   1  'Fixed Single
         Height          =   855
         Left            =   1020
         TabIndex        =   91
         Top             =   690
         Width           =   885
      End
      Begin VB.Label Label71 
         BackColor       =   &H00000000&
         Caption         =   "Right Resp"
         ForeColor       =   &H0000FFFF&
         Height          =   465
         Left            =   2160
         TabIndex        =   68
         Top             =   210
         Width           =   495
      End
      Begin VB.Label Label70 
         BackColor       =   &H00000000&
         Caption         =   "Left Resp"
         ForeColor       =   &H0000FFFF&
         Height          =   465
         Left            =   1260
         TabIndex        =   67
         Top             =   210
         Width           =   495
      End
      Begin VB.Label Label69 
         BackColor       =   &H00000000&
         Caption         =   "Right Trial"
         ForeColor       =   &H0000FFFF&
         Height          =   465
         Left            =   420
         TabIndex        =   66
         Top             =   1830
         Width           =   495
      End
      Begin VB.Label Label68 
         BackColor       =   &H00000000&
         Caption         =   "Left Trial"
         ForeColor       =   &H0000FFFF&
         Height          =   465
         Left            =   420
         TabIndex        =   65
         Top             =   930
         Width           =   495
      End
   End
   Begin Mh3dfrmLibCtl.Mh3dFrame Mh3dFrame2 
      Height          =   3225
      Left            =   3360
      TabIndex        =   69
      Top             =   600
      Width           =   3495
      _Version        =   65536
      _ExtentX        =   6165
      _ExtentY        =   5689
      _StockProps     =   79
      Caption         =   "Mh3dFrame2"
      BeginProperty Font {0BE35203-8F91-11CE-9DE3-00AA004BB851} 
         Name            =   "MS Sans Serif"
         Size            =   8.25
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Alignment       =   0
      AutoSize        =   0   'False
      BevelSize       =   4
      BevelStyle      =   0
      BorderColor     =   -2147483642
      BorderStyle     =   0
      FillColor       =   0
      FontStyle       =   0
      FontTransparent =   0   'False
      LightColor      =   -2147483643
      ShadowColor     =   -2147483632
      TextColor       =   -2147483640
      WallPaper       =   0
      NoPrefix        =   0   'False
      FormatString    =   ""
      Caption         =   "Mh3dFrame2"
      Picture         =   "CONTROL.frx":B83C
      Begin LeadLib.Lead leadStimuli 
         Height          =   1095
         Left            =   240
         TabIndex        =   79
         Top             =   240
         Width           =   1335
         _Version        =   393216
         _ExtentX        =   2355
         _ExtentY        =   1931
         _StockProps     =   97
         BorderStyle     =   1
         ScaleHeight     =   1065
         ScaleMode       =   1
         ScaleWidth      =   1305
         DataField       =   ""
      End
      Begin LeadLib.Lead leadProbe 
         Height          =   1095
         Left            =   240
         TabIndex        =   80
         Top             =   1800
         Width           =   1335
         _Version        =   393216
         _ExtentX        =   2355
         _ExtentY        =   1931
         _StockProps     =   97
         BorderStyle     =   1
         ScaleHeight     =   1065
         ScaleMode       =   1
         ScaleWidth      =   1305
         DataField       =   ""
      End
      Begin LeadLib.Lead leadProbe2 
         Height          =   1095
         Left            =   1920
         TabIndex        =   114
         Top             =   1080
         Width           =   1335
         _Version        =   393216
         _ExtentX        =   2355
         _ExtentY        =   1931
         _StockProps     =   97
         BorderStyle     =   1
         ScaleHeight     =   1065
         ScaleMode       =   1
         ScaleWidth      =   1305
         DataField       =   ""
      End
      Begin VB.Label lblDiffSelection 
         Alignment       =   1  'Right Justify
         BackColor       =   &H00FFFFFF&
         Height          =   495
         Left            =   2880
         TabIndex        =   70
         Top             =   240
         Visible         =   0   'False
         Width           =   525
      End
   End
   Begin VB.Label Label52 
      BackColor       =   &H00404040&
      Caption         =   "Probe Abort LTO:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   6540
      TabIndex        =   113
      Top             =   7380
      Visible         =   0   'False
      Width           =   1245
   End
   Begin VB.Label Label51 
      BackColor       =   &H00404040&
      Caption         =   "msec."
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   9000
      TabIndex        =   112
      Top             =   7350
      Visible         =   0   'False
      Width           =   405
   End
   Begin VB.Label Label50 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "ProbeAbortLTO"
      DataSource      =   "datSessionInfo"
      Height          =   285
      Left            =   8010
      TabIndex        =   111
      Top             =   7320
      Visible         =   0   'False
      Width           =   885
   End
   Begin VB.Label Label49 
      BackColor       =   &H00404040&
      Caption         =   "Probe Aborts Enabled:"
      ForeColor       =   &H00C0C0C0&
      Height          =   225
      Left            =   6540
      TabIndex        =   110
      Top             =   7080
      Visible         =   0   'False
      Width           =   2205
   End
   Begin VB.Label Label46 
      BackColor       =   &H00404040&
      Caption         =   "Touch Aborts Enabled:"
      ForeColor       =   &H00C0C0C0&
      Height          =   225
      Left            =   3630
      TabIndex        =   108
      Top             =   6810
      Visible         =   0   'False
      Width           =   2205
   End
   Begin VB.Label Label45 
      BackColor       =   &H00404040&
      Caption         =   "FR - 0   Session:"
      ForeColor       =   &H00C0C0C0&
      Height          =   225
      Left            =   3780
      TabIndex        =   106
      Top             =   6720
      Visible         =   0   'False
      Width           =   2205
   End
   Begin VB.Label lblTransfer 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "PercentTCorrect"
      DataSource      =   "datSessionResults"
      Height          =   255
      Left            =   8520
      TabIndex        =   104
      Top             =   5400
      Width           =   885
   End
   Begin VB.Label Label44 
      BackColor       =   &H00404040&
      Caption         =   "% Transfer:"
      ForeColor       =   &H00C0C0C0&
      Height          =   225
      Left            =   7260
      TabIndex        =   103
      Top             =   5430
      Width           =   1155
   End
   Begin VB.Label lblSetSize 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "SetSize"
      DataSource      =   "datSessionInfo"
      Height          =   285
      Left            =   1650
      TabIndex        =   102
      Top             =   5340
      Width           =   885
   End
   Begin VB.Label Label43 
      BackColor       =   &H00404040&
      Caption         =   "Set Size:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   180
      TabIndex        =   101
      Top             =   5340
      Width           =   1305
   End
   Begin VB.Label Label42 
      BackColor       =   &H00404040&
      Caption         =   "Fixed Ratio (FR):"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   180
      TabIndex        =   100
      Top             =   4950
      Width           =   1305
   End
   Begin VB.Label Label37 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "StimTouches"
      DataSource      =   "datSessionInfo"
      Height          =   285
      Left            =   1650
      TabIndex        =   99
      Top             =   4920
      Width           =   885
   End
   Begin VB.Label lblMode 
      Alignment       =   2  'Center
      BackColor       =   &H0000C000&
      BorderStyle     =   1  'Fixed Single
      Caption         =   "Normal Mode"
      BeginProperty Font 
         Name            =   "MS Sans Serif"
         Size            =   8.25
         Charset         =   0
         Weight          =   700
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      ForeColor       =   &H00000000&
      Height          =   315
      Left            =   3360
      TabIndex        =   90
      Top             =   0
      Width           =   3495
   End
   Begin VB.Label lblResponse 
      Alignment       =   2  'Center
      BackColor       =   &H00000000&
      BorderStyle     =   1  'Fixed Single
      Caption         =   "Same/SameTrial"
      BeginProperty Font 
         Name            =   "MS Sans Serif"
         Size            =   12
         Charset         =   0
         Weight          =   700
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      ForeColor       =   &H0000FFFF&
      Height          =   405
      Left            =   7290
      TabIndex        =   89
      Top             =   1530
      Width           =   2145
   End
   Begin VB.Label lblApparatus 
      Alignment       =   2  'Center
      BackColor       =   &H00000000&
      BorderStyle     =   1  'Fixed Single
      Caption         =   "Grain"
      BeginProperty Font 
         Name            =   "MS Sans Serif"
         Size            =   18
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      ForeColor       =   &H0000FFFF&
      Height          =   525
      Left            =   7230
      TabIndex        =   88
      Top             =   480
      Width           =   2145
   End
   Begin VB.Label Label38 
      BackColor       =   &H00404040&
      Caption         =   "Juice Delay:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   7260
      TabIndex        =   87
      Top             =   2280
      Visible         =   0   'False
      Width           =   945
   End
   Begin VB.Label Label39 
      BackColor       =   &H00404040&
      Caption         =   "Grain Delay:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   7260
      TabIndex        =   86
      Top             =   2610
      Width           =   945
   End
   Begin VB.Label Label40 
      BackColor       =   &H00404040&
      Caption         =   "msec."
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   8910
      TabIndex        =   84
      Top             =   2280
      Visible         =   0   'False
      Width           =   525
   End
   Begin VB.Label Label41 
      BackColor       =   &H00404040&
      Caption         =   "msec."
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   8910
      TabIndex        =   82
      Top             =   2610
      Width           =   525
   End
   Begin VB.Label lblCorrectPosition 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "Configuration"
      DataSource      =   "datTrialInfo"
      Height          =   285
      Left            =   8520
      TabIndex        =   78
      Top             =   4500
      Width           =   885
   End
   Begin VB.Label Label14 
      BackColor       =   &H00404040&
      Caption         =   "Configuration (DS):"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   7020
      TabIndex        =   77
      Top             =   4530
      Width           =   1395
   End
   Begin VB.Label lblTotalIncorrectCorrections 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      Height          =   285
      Left            =   8520
      TabIndex        =   76
      Top             =   6180
      Width           =   885
   End
   Begin VB.Label lblTotalAbortCorrections 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      Height          =   285
      Left            =   9240
      TabIndex        =   75
      Top             =   7050
      Visible         =   0   'False
      Width           =   885
   End
   Begin VB.Label Label35 
      BackColor       =   &H00404040&
      Caption         =   "Abort Corr.:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   7980
      TabIndex        =   74
      Top             =   7080
      Visible         =   0   'False
      Width           =   1155
   End
   Begin VB.Label lblNumIncorrectCorrections 
      BackColor       =   &H00404040&
      Caption         =   "Incorrect Corr.:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   7260
      TabIndex        =   73
      Top             =   6210
      Width           =   1155
   End
   Begin VB.Label lblClock 
      Alignment       =   2  'Center
      BorderStyle     =   1  'Fixed Single
      Caption         =   "00:00:00"
      Height          =   285
      Left            =   2550
      TabIndex        =   72
      Top             =   390
      Width           =   765
   End
   Begin VB.Label lblStateMessage 
      BorderStyle     =   1  'Fixed Single
      Height          =   285
      Left            =   3360
      TabIndex        =   71
      Top             =   300
      Width           =   3495
   End
   Begin VB.Label Label67 
      BackColor       =   &H00404040&
      Caption         =   "Response"
      BeginProperty Font 
         Name            =   "MS Sans Serif"
         Size            =   13.5
         Charset         =   0
         Weight          =   700
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      ForeColor       =   &H00C0C0C0&
      Height          =   375
      Left            =   7230
      TabIndex        =   63
      Top             =   1140
      Width           =   1515
   End
   Begin VB.Label Label66 
      BackColor       =   &H00404040&
      Caption         =   "Apparatus"
      BeginProperty Font 
         Name            =   "MS Sans Serif"
         Size            =   13.5
         Charset         =   0
         Weight          =   700
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      ForeColor       =   &H00C0C0C0&
      Height          =   375
      Left            =   7230
      TabIndex        =   62
      Top             =   90
      Width           =   1515
   End
   Begin VB.Label Label65 
      BackColor       =   &H00404040&
      Caption         =   "Summary"
      BeginProperty Font 
         Name            =   "MS Sans Serif"
         Size            =   13.5
         Charset         =   0
         Weight          =   700
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      ForeColor       =   &H00C0C0C0&
      Height          =   345
      Left            =   7620
      TabIndex        =   61
      Top             =   3000
      Width           =   1515
   End
   Begin VB.Label Label63 
      BackColor       =   &H00404040&
      Caption         =   "# Grains:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   7260
      TabIndex        =   60
      Top             =   6630
      Width           =   735
   End
   Begin VB.Label Label62 
      BackColor       =   &H00404040&
      Caption         =   "# Timeouts:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   7260
      TabIndex        =   59
      Top             =   5820
      Width           =   855
   End
   Begin VB.Label Label61 
      BackColor       =   &H00404040&
      Caption         =   "# Aborts:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   8100
      TabIndex        =   58
      Top             =   6840
      Visible         =   0   'False
      Width           =   795
   End
   Begin VB.Label Label60 
      BackColor       =   &H00404040&
      Caption         =   "% Correct:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   7260
      TabIndex        =   57
      Top             =   4980
      Width           =   825
   End
   Begin VB.Label Label59 
      BackColor       =   &H00404040&
      Caption         =   "Trial #:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   7260
      TabIndex        =   56
      Top             =   3630
      Width           =   615
   End
   Begin VB.Label Label58 
      BackColor       =   &H00404040&
      Caption         =   "Type:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   7260
      TabIndex        =   55
      Top             =   4080
      Width           =   855
   End
   Begin VB.Label lblType 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "TrialType"
      DataSource      =   "datTrialInfo"
      Height          =   285
      Left            =   8520
      TabIndex        =   54
      Top             =   4050
      Width           =   885
   End
   Begin VB.Label lblPercentCorrect 
      Alignment       =   1  'Right Justify
      BackColor       =   &H000000FF&
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "PercentCorrect"
      DataSource      =   "datSessionResults"
      Height          =   285
      Left            =   8520
      TabIndex        =   53
      Top             =   4950
      Width           =   885
   End
   Begin VB.Label lblTrialNum 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "TrialNum"
      DataSource      =   "datTrialInfo"
      Height          =   285
      Left            =   8520
      TabIndex        =   52
      Top             =   3600
      Width           =   885
   End
   Begin VB.Label lblNumAborts 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      Height          =   285
      Left            =   9360
      TabIndex        =   51
      Top             =   6780
      Visible         =   0   'False
      Width           =   885
   End
   Begin VB.Label lblNumJuice 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      Height          =   285
      Left            =   8520
      TabIndex        =   50
      Top             =   5790
      Width           =   885
   End
   Begin VB.Label lblNumPellets 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      Height          =   285
      Left            =   8520
      TabIndex        =   49
      Top             =   6600
      Width           =   885
   End
   Begin VB.Label lblDate 
      Alignment       =   2  'Center
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      Height          =   285
      Left            =   810
      TabIndex        =   48
      Top             =   390
      Width           =   1725
   End
   Begin VB.Label lblIncorrectLTO 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "DefaultIncorrectLTO"
      DataSource      =   "datSessionInfo"
      Height          =   285
      Left            =   1650
      TabIndex        =   47
      Top             =   6090
      Width           =   885
   End
   Begin VB.Label lblAbortLTO 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "DefaultAbortLTO"
      DataSource      =   "datSessionInfo"
      Height          =   285
      Left            =   5100
      TabIndex        =   46
      Top             =   7020
      Visible         =   0   'False
      Width           =   885
   End
   Begin VB.Label lblDiffProbability 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "DiffProbability"
      DataSource      =   "datSessionInfo"
      Height          =   285
      Left            =   1650
      TabIndex        =   45
      Top             =   4500
      Width           =   885
   End
   Begin VB.Label lblSameProbability 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "SameProbability"
      DataSource      =   "datSessionInfo"
      Height          =   285
      Left            =   1650
      TabIndex        =   44
      Top             =   4110
      Width           =   885
   End
   Begin VB.Label lblReinfProbability 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "ReinfProbability"
      DataSource      =   "datSessionInfo"
      Height          =   285
      Left            =   1650
      TabIndex        =   43
      Top             =   3720
      Width           =   885
   End
   Begin VB.Label lblPelletProbability 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "PelletProbability"
      DataSource      =   "datSessionInfo"
      Height          =   285
      Left            =   1650
      TabIndex        =   42
      Top             =   3330
      Width           =   885
   End
   Begin VB.Label lblSessionID 
      Alignment       =   2  'Center
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "SessionID"
      DataSource      =   "datSessionInfo"
      Height          =   285
      Left            =   1650
      TabIndex        =   41
      Top             =   1320
      Width           =   885
   End
   Begin VB.Label lblSubject 
      Alignment       =   2  'Center
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "Subject"
      DataSource      =   "datSessionInfo"
      Height          =   285
      Left            =   1650
      TabIndex        =   40
      Top             =   870
      Width           =   885
   End
   Begin VB.Label lblNumOfTrials 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "NumOfTrials"
      DataSource      =   "datSessionInfo"
      Height          =   285
      Left            =   1650
      TabIndex        =   39
      Top             =   1770
      Width           =   885
   End
   Begin VB.Label lblListLength 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "ListLength"
      DataSource      =   "datTrialInfo"
      Height          =   285
      Left            =   4890
      TabIndex        =   38
      Top             =   6270
      Visible         =   0   'False
      Width           =   885
   End
   Begin VB.Label lblViewTime 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      Height          =   285
      Left            =   4890
      TabIndex        =   37
      Top             =   6660
      Visible         =   0   'False
      Width           =   885
   End
   Begin VB.Label lblITI 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "ITI"
      DataSource      =   "datTrialInfo"
      Height          =   285
      Left            =   1650
      TabIndex        =   36
      Top             =   2130
      Width           =   885
   End
   Begin VB.Label lblISI 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      Height          =   285
      Left            =   1650
      TabIndex        =   35
      Top             =   2520
      Width           =   885
   End
   Begin VB.Label lblJuiceProbability 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "JuiceProbability"
      DataSource      =   "datSessionInfo"
      Height          =   285
      Left            =   5130
      TabIndex        =   34
      Top             =   7260
      Visible         =   0   'False
      Width           =   885
   End
   Begin VB.Label lblProbeDelay 
      Alignment       =   1  'Right Justify
      BorderStyle     =   1  'Fixed Single
      Caption         =   "0"
      DataField       =   "Delay"
      DataSource      =   "datProbeInfo"
      Height          =   285
      Left            =   1650
      TabIndex        =   33
      Top             =   2910
      Width           =   885
   End
   Begin VB.Label Label1 
      BackColor       =   &H00404040&
      Caption         =   "Date:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   180
      TabIndex        =   32
      Top             =   420
      Width           =   525
   End
   Begin VB.Label Label16 
      BackColor       =   &H00404040&
      Caption         =   "sec."
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   2640
      TabIndex        =   31
      Top             =   6150
      Width           =   405
   End
   Begin VB.Label Label15 
      BackColor       =   &H00404040&
      Caption         =   "sec."
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   6090
      TabIndex        =   30
      Top             =   7050
      Visible         =   0   'False
      Width           =   405
   End
   Begin VB.Label Label13 
      BackColor       =   &H00404040&
      Caption         =   "Choice Delay:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   180
      TabIndex        =   29
      Top             =   2940
      Width           =   1005
   End
   Begin VB.Label Label12 
      BackColor       =   &H00404040&
      Caption         =   "ISI:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   180
      TabIndex        =   28
      Top             =   2520
      Width           =   345
   End
   Begin VB.Label Label11 
      BackColor       =   &H00404040&
      Caption         =   "ITI:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   180
      TabIndex        =   27
      Top             =   2130
      Width           =   345
   End
   Begin VB.Label Label10 
      BackColor       =   &H00404040&
      Caption         =   "View Time:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   3420
      TabIndex        =   26
      Top             =   6660
      Visible         =   0   'False
      Width           =   885
   End
   Begin VB.Label Label9 
      BackColor       =   &H00404040&
      Caption         =   "List Length:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   3420
      TabIndex        =   25
      Top             =   6270
      Visible         =   0   'False
      Width           =   915
   End
   Begin VB.Label Label8 
      BackColor       =   &H00404040&
      Caption         =   "Time Out:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   180
      TabIndex        =   24
      Top             =   6120
      Width           =   1065
   End
   Begin VB.Label Label7 
      BackColor       =   &H00404040&
      Caption         =   "Abort LTO:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   3630
      TabIndex        =   23
      Top             =   7080
      Visible         =   0   'False
      Width           =   1065
   End
   Begin VB.Label Label6 
      BackColor       =   &H00404040&
      Caption         =   "Correction Procedure:"
      ForeColor       =   &H00C0C0C0&
      Height          =   225
      Left            =   180
      TabIndex        =   22
      Top             =   5730
      Width           =   2205
   End
   Begin VB.Label Label5 
      BackColor       =   &H00404040&
      Caption         =   "Abort Corrections Enabled:"
      ForeColor       =   &H00C0C0C0&
      Height          =   225
      Left            =   3660
      TabIndex        =   21
      Top             =   6930
      Visible         =   0   'False
      Width           =   2055
   End
   Begin VB.Label Label4 
      BackColor       =   &H00404040&
      Caption         =   "Number of Trials:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   180
      TabIndex        =   20
      Top             =   1800
      Width           =   1275
   End
   Begin VB.Label Label3 
      BackColor       =   &H00404040&
      Caption         =   "Session ID:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   180
      TabIndex        =   19
      Top             =   1350
      Width           =   915
   End
   Begin VB.Label Label2 
      BackColor       =   &H00404040&
      Caption         =   "Subject:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   180
      TabIndex        =   18
      Top             =   900
      Width           =   675
   End
   Begin VB.Label Label17 
      BackColor       =   &H00404040&
      Caption         =   "sec."
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   2640
      TabIndex        =   17
      Top             =   2940
      Width           =   405
   End
   Begin VB.Label Label18 
      BackColor       =   &H00404040&
      Caption         =   "sec."
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   2640
      TabIndex        =   16
      Top             =   2550
      Width           =   405
   End
   Begin VB.Label Label19 
      BackColor       =   &H00404040&
      Caption         =   "sec."
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   2640
      TabIndex        =   15
      Top             =   2130
      Width           =   405
   End
   Begin VB.Label Label20 
      BackColor       =   &H00404040&
      Caption         =   "msec."
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   5880
      TabIndex        =   14
      Top             =   6660
      Visible         =   0   'False
      Width           =   405
   End
   Begin VB.Label Label21 
      BackColor       =   &H00404040&
      Caption         =   "Juice Probability:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   3660
      TabIndex        =   13
      Top             =   7290
      Visible         =   0   'False
      Width           =   1305
   End
   Begin VB.Label Label22 
      BackColor       =   &H00404040&
      Caption         =   "Grain Probability:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   180
      TabIndex        =   12
      Top             =   3360
      Width           =   1305
   End
   Begin VB.Label Label23 
      BackColor       =   &H00404040&
      Caption         =   "Reinf Probability:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   180
      TabIndex        =   11
      Top             =   3750
      Width           =   1305
   End
   Begin VB.Label Label24 
      BackColor       =   &H00404040&
      Caption         =   "Left Probability:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   180
      TabIndex        =   10
      Top             =   4140
      Width           =   1305
   End
   Begin VB.Label Label25 
      BackColor       =   &H00404040&
      Caption         =   "Right Probability:"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   180
      TabIndex        =   9
      Top             =   4530
      Width           =   1305
   End
   Begin VB.Label Label26 
      BackColor       =   &H00404040&
      Caption         =   "%"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   6120
      TabIndex        =   8
      Top             =   7290
      Visible         =   0   'False
      Width           =   345
   End
   Begin VB.Label Label27 
      BackColor       =   &H00404040&
      Caption         =   "%"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   2640
      TabIndex        =   7
      Top             =   3360
      Width           =   345
   End
   Begin VB.Label Label28 
      BackColor       =   &H00404040&
      Caption         =   "%"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   2640
      TabIndex        =   6
      Top             =   3750
      Width           =   345
   End
   Begin VB.Label Label29 
      BackColor       =   &H00404040&
      Caption         =   "%"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   2640
      TabIndex        =   5
      Top             =   4140
      Width           =   345
   End
   Begin VB.Label Label30 
      BackColor       =   &H00404040&
      Caption         =   "%"
      ForeColor       =   &H00C0C0C0&
      Height          =   255
      Left            =   2640
      TabIndex        =   4
      Top             =   4530
      Width           =   345
   End
   Begin VB.Label Label31 
      BackColor       =   &H00404040&
      Caption         =   "Variable Probe Delay"
      ForeColor       =   &H00C0C0C0&
      Height          =   225
      Left            =   6540
      TabIndex        =   3
      Top             =   6870
      Visible         =   0   'False
      Width           =   1515
   End
   Begin VB.Menu mnuFile 
      Caption         =   "&File"
      Begin VB.Menu mnuClose 
         Caption         =   "&Close"
      End
   End
   Begin VB.Menu mnuSpace2 
      Caption         =   ""
   End
   Begin VB.Menu mnuNormalMode 
      Caption         =   "&Normal"
   End
   Begin VB.Menu mnuTestMode 
      Caption         =   "&Test"
   End
   Begin VB.Menu mnuSpace1 
      Caption         =   ""
   End
   Begin VB.Menu mnuStart 
      Caption         =   "&Start"
   End
   Begin VB.Menu mnuPause 
      Caption         =   "&Pause"
      Enabled         =   0   'False
   End
   Begin VB.Menu mnuResume 
      Caption         =   "&Resume"
      Enabled         =   0   'False
   End
   Begin VB.Menu mnuReset 
      Caption         =   "R&eset"
      Enabled         =   0   'False
   End
End
Attribute VB_Name = "frmControl"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False
Option Explicit
 
Dim sFile, sFile1, sFile2 As String
Dim Pos1, Pos2 As Integer

Const PROBE_LOCATION = 13

'loading stim files as an array will simply things later on.
'stimFile PROBE_LOCATION is ALWAYS the "change-to" image
'stimFile 1 is ALWAYS the "changed" item
'stimFile 2 is ALWAYS the "no change" item
Dim stimulusFiles(13) As String

'storing the positions in an array facilitates looping later on
Dim positions(12) As Integer


Dim lResponseTime As Long
Dim sResponse As String
Dim RespLoc As Integer
Dim HoldResponse As String
Dim sTimerState As String
Dim sInitState As String
Dim sImageDisplayed As String
Dim iStimTouches As Integer
Dim iProbeTouches As Integer
Dim gPecks, gDisplayTime As Long
Dim Item1 As String
Dim Item2 As String
Dim Item3 As String
Dim Item4 As String
Dim Item5 As String
Dim Item6 As String
Dim Item7 As String
Dim Item8 As String
Dim Item9 As String
Dim Item10 As String
Dim ItemP As String
Dim ItemP2 As String
Dim butt As Integer
Dim flagAbort As Boolean
Dim TrialStatus As Boolean


Public Sub showReadySignal()

    With frmDisplay.leadProbe
        .Load "C:/monkey_cdt_images_1/rs.jpg", 0, 0
        .ForceRepaint
        .DstWidth = .Width
        .DstHeight = .Height
        .Visible = True
    End With

'    Debug.Print " Showing RS"

End Sub


Public Sub handleReadySignalTouch()
    frmDisplay.leadProbe.Visible = False
    
    tmrAfterReady.Interval = 500
    tmrAfterReady.Enabled = True
End Sub

Private Sub tmrAfterReady_Timer()
    tmrAfterReady.Enabled = False
    subInitialItems
End Sub







Public Sub DetectedRightPeck()
    'sResponse = "Right"
    'subCheckResponse
End Sub

Public Sub DetectedFormPeck()
    If listenToFormPecks = True Then
        RespLoc = 0
        sResponse = "Form"
        'If sImageDisplayed = "first" Then subRecordTouch
        'subCheckResponse
    End If
End Sub

Public Sub DetectedLeftPeck()
    'sResponse = "Left"
    'subCheckResponse
End Sub

Public Sub DetectedPos1()
        RespLoc = 1
        If sImageDisplayed = "first" Then subRecordTouch
        subCheckResponse
End Sub
Public Sub DetectedPos2()
        RespLoc = 2
        If sImageDisplayed = "first" Then subRecordTouch
        subCheckResponse
End Sub
Public Sub DetectedPos3()
        RespLoc = 3
        If sImageDisplayed = "first" Then subRecordTouch
        subCheckResponse
End Sub
Public Sub DetectedPos4()
        RespLoc = 4
        If sImageDisplayed = "first" Then subRecordTouch
        subCheckResponse
End Sub
Public Sub DetectedPos5()
        RespLoc = 5
        If sImageDisplayed = "first" Then subRecordTouch
        subCheckResponse
End Sub
Public Sub DetectedPos6()
        RespLoc = 6
        If sImageDisplayed = "first" Then subRecordTouch
        subCheckResponse
End Sub
Public Sub DetectedPos7()
        RespLoc = 7
        If sImageDisplayed = "first" Then subRecordTouch
        subCheckResponse
End Sub
Public Sub DetectedPos8()
        RespLoc = 8
        If sImageDisplayed = "first" Then subRecordTouch
        subCheckResponse
End Sub
Public Sub DetectedPos9()
        RespLoc = 9
        If sImageDisplayed = "first" Then subRecordTouch
        subCheckResponse
End Sub
Public Sub DetectedPos10()
        RespLoc = 10
        If sImageDisplayed = "first" Then subRecordTouch
        subCheckResponse
End Sub
Public Sub DetectedPos11()
        RespLoc = 11
        If sImageDisplayed = "first" Then subRecordTouch
        subCheckResponse
End Sub
Public Sub DetectedPos12()
        RespLoc = 12
        If sImageDisplayed = "first" Then subRecordTouch
        subCheckResponse
End Sub
Public Sub DetectedPos13()
        RespLoc = 13
        If sImageDisplayed = "first" Then subRecordTouch
        subCheckResponse
End Sub
Public Sub DetectedPos14()
        RespLoc = 14
        If sImageDisplayed = "first" Then subRecordTouch
        subCheckResponse
End Sub
Public Sub DetectedPos15()
        RespLoc = 15
        If sImageDisplayed = "first" Then subRecordTouch
        subCheckResponse
End Sub
Public Sub DetectedPos16()
        RespLoc = 16
        If sImageDisplayed = "first" Then subRecordTouch
        subCheckResponse
End Sub

Public Sub DetectedSamplePeck()
        'subRecordTouch
        'If sImageDisplayed = "stim" Then sResponse = "Stim"
        'subCheckResponse
End Sub
Private Sub subTranPercentCorrect()
                                                    
    Dim iDenominator As Integer
    
    datSessionResults.Recordset.MoveFirst
    datSessionResults.Recordset.Edit
    iDenominator = datSessionResults.Recordset("TSameCorrect") _
            + datSessionResults.Recordset("TSameIncorrect") _
            + datSessionResults.Recordset("TDiffCorrect") _
            + datSessionResults.Recordset("TDiffIncorrect")
       
    If iDenominator = 0 Then
        iDenominator = 1
    End If
    
    datSessionResults.Recordset("PercentTCorrect") = Format(100 * (datSessionResults.Recordset("TSameCorrect") _
                                + datSessionResults.Recordset("TDiffCorrect")) _
                                / iDenominator, "##0")
    datSessionResults.Recordset.Update

End Sub
Private Sub subCheckResponse()
'    cursor_visible (True)
    Dim FR As Integer
    FR = datSessionInfo.Recordset("StimTouches")
    Dim sTrialType As String
    Dim iTotalStimTouches As Integer
    Dim lISITime As Long
    Dim NoPeckTime As Integer
    
    'DRO fixed at 1s
    NoPeckTime = 1000 'milliseconds
    
    If sImageDisplayed = "first" Then
        
        iTotalStimTouches = datSessionInfo.Recordset("StimTouches")
        
        Dim useDRO As Boolean
        If datTrialInfo.Recordset("TrialType") = "DRO" Then
            useDRO = True
        Else
            useDRO = False
        End If
        
        If sResponse = "Stim" Then
            Stop
            iStimTouches = iStimTouches + 1
            lblStateMessage.Caption = "Stim Touches = " & iStimTouches & " / " & iTotalStimTouches
            sResponse = "None"
        End If
        
        If (lResponseTime >= gDisplayTime) And (gPecks >= FR) Then
            
            sImageDisplayed = "idle"
            tmrProbeView.Enabled = False
            tmrProbeView.Interval = 0
            subClearDisplay
            lISITime = datStimuliInfo.Recordset("ISI")
            lblStateMessage.Caption = "Probe Delay (" & datProbeInfo.Recordset("Delay") & " sec.)"
            
            If (datProbeInfo.Recordset("Delay") = 0) Then
                subSecondItems
            Else
                tmrProbeDelay.Interval = datProbeInfo.Recordset("Delay")
                tmrProbeDelay.Enabled = True
            End If
            
        End If
        
        Debug.Print ("pecks: " & gPecks & ", time: " & lResponseTime & ", vt: " & gDisplayTime)
'        If (lResponseTime >= (gDisplayTime * 1000)) And (gPecks > 0) Then
        If (useDRO And gPecks > 0) Then
            gDisplayTime = gDisplayTime + NoPeckTime
            gPecks = 0
        End If
                                            'choice after 500 ms
    ElseIf (sImageDisplayed = "second") And RespLoc <> 0 Then ' (lResponseTime >= 500) Then
                    
        If datTrialInfo.Recordset("POSchange") = RespLoc Then 'correct
            
            listenToFormPecks = False
        
            tmrProbeView.Enabled = False
            tmrProbeView.Interval = 0
            subClearDisplay
            If lblTotalIncorrectCorrections.Caption = 0 Then lblSameTrialSameResp.Caption = lblSameTrialSameResp.Caption + 1
            TrialStatus = True
            subCorrect
        End If
        If datTrialInfo.Recordset("POSchange") <> RespLoc Then  'incorrect
            
            listenToFormPecks = False
            
            tmrProbeView.Enabled = False
            tmrProbeView.Interval = 0
            subClearDisplay
            If lblTotalIncorrectCorrections.Caption = 0 Then lblSameTrialDiffResp.Caption = lblSameTrialDiffResp.Caption + 1
            TrialStatus = False
            If datSessionInfo.Recordset("IncorrectCorrectionsEnabled") = True Then
                subIncorrectLTO
            Else
                subIncorrect
            End If
        End If
    
    End If
'    cursor_visible (False)
End Sub
Private Sub subClearDisplay()
    

        frmDisplay.leadStimuli.Visible = False
        frmDisplay.Lead1.Visible = False
        frmDisplay.Lead2.Visible = False
        frmDisplay.Lead3.Visible = False
        frmDisplay.Lead4.Visible = False
        frmDisplay.Lead5.Visible = False
        frmDisplay.Lead6.Visible = False
        frmDisplay.Lead7.Visible = False
        frmDisplay.Lead8.Visible = False
        frmDisplay.Lead9.Visible = False
        frmDisplay.Lead10.Visible = False
        frmDisplay.Lead11.Visible = False
        frmDisplay.Lead12.Visible = False
        frmDisplay.Lead13.Visible = False
        frmDisplay.Lead14.Visible = False
        frmDisplay.Lead15.Visible = False
        leadStimuli.Visible = False
        leadProbe.Visible = False
        leadProbe2.Visible = False
        
  '  If TrialStatus = True Then
        'agOutp METRABYTE_CONTROL, APPARATUS_SET_B_OUTPUT
        'agOutp METRABYTE_PORT_B, APPARATUS_HOUSE_LIGHT
'With DriverLINXSR1
'    .Req_subsystem = DL_DO
'    .Req_mode = DL_POLLED
'    .Req_op = DL_START
'    .Evt_Str_type = DL_NULLEVENT
'    .Evt_Stp_type = DL_NULLEVENT
'    .Evt_Tim_type = DL_NULLEVENT
'    .Sel_chan_format = DL_tNATIVE
'    .Sel_chan_N = 1
'    .Sel_chan_start = 1   'Channel 0 is configured for Output already
'    .Sel_chan_startGainCode = 0
'    .Sel_buf_samples = 1
'    .Sel_buf_N = 0
'    .Refresh
'    .Res_Sta_ioValue = APPARATUS_HOUSE_LIGHT
'    .Refresh
'    Label1.Caption = .Message
'End With
        
        
 '       lblStateMessage.Caption = "ITI (" & datTrialInfo.Recordset("ITI") & " sec.)"
 '       tmrITI.Interval = datTrialInfo.Recordset("ITI") * 1000
 '       flagAbort = False
 '       tmrITI.Enabled = True
 '   End If
    
 '   If TrialStatus = False Then
 '       lblStateMessage.Caption = "Time Out (" & datSessionInfo.Recordset("DefaultIncorrectLTO") & " sec.)"
'        lblTotalIncorrectCorrections.Caption = lblTotalIncorrectCorrections.Caption + 1
 '       tmrIncorrectLTO.Interval = datSessionInfo.Recordset("DefaultIncorrectLTO") * 1000
  '      If tmrIncorrectLTO.Interval = 0 Then tmrIncorrectLTO.Interval = 1
  '      tmrIncorrectLTO.Enabled = True
  '  End If

End Sub

Private Sub subDisplayTraining()

    Dim sFile As String
  '  Dim lViewTime As Long
  '
  '  On Error GoTo Error_subDisplayTraining'''

    
   ' lResponseTime = 0
   ' sResponse = "None"
   ' lblResponse.Caption = ""
   ' lblApparatus.Caption = ""
   ' lblStateMessage.Caption = "Training"
    
    'Stimuli
  '  sFile = datSessionInfo.Recordset("ImagePath") & "\" & datStimuliInfo.Recordset("Filename")
  '
  '  leadStimuli.AutoScroll = False
  '  leadStimuli.Load sFile, 0, 0
  '  leadStimuli.DstTop = 0
  '  leadStimuli.DstLeft = 0
  '  leadStimuli.DstWidth = leadStimuli.Width
  '  leadStimuli.DstHeight = leadStimuli.Height
  '  leadStimuli.ForceRepaint
  '  leadStimuli.Visible = True
  '
   ' If datSessionInfo.Recordset("OperationalMode") = "Normal" Then
  '      frmDisplay.leadStimuli.Load sFile, 0, 0
  '      frmDisplay.leadStimuli.DstTop = 0
  '      frmDisplay.leadStimuli.DstLeft = 0
  '      frmDisplay.leadStimuli.DstWidth = frmDisplay.leadStimuli.Width
  '      frmDisplay.leadStimuli.DstHeight = frmDisplay.leadStimuli.Height
  '      frmDisplay.leadStimuli.ForceRepaint
  '      frmDisplay.leadStimuli.Visible = True
  '  End If
    
    
    'Probe
  'mt  gsubEnableTouchscreen tmrReadTouch
    
   ' sFile = datSessionInfo.Recordset("ImagePath") & "\" & datProbeInfo.Recordset("Filename")

   ' leadProbe.AutoScroll = False
   ' leadProbe.Load sFile, 0, 0
   ' leadProbe.DstTop = 0
   ' leadProbe.DstLeft = 0
   ' leadProbe.DstWidth = leadProbe.Width
   ' leadProbe.DstHeight = leadProbe.Height
   ' leadProbe.ForceRepaint
   ' leadProbe.Visible = True
   '
   ' If datSessionInfo.Recordset("OperationalMode") = "Test" Then
       ' lblDiffSelection.Visible = True
        
  '  ElseIf datSessionInfo.Recordset("OperationalMode") = "Normal" Then
  '      frmDisplay.leadProbe.Load sFile, 0, 0
  '      frmDisplay.leadProbe.DstTop = 0
  '      frmDisplay.leadProbe.DstLeft = 0
  '      frmDisplay.leadProbe.DstWidth = frmDisplay.leadProbe.Width
   '     frmDisplay.leadProbe.DstHeight = frmDisplay.leadProbe.Height
   '     frmDisplay.leadProbe.ForceRepaint
   '     frmDisplay.leadProbe.Visible = True
        'frmDisplay.lblDiff.Visible = True
    
  '      frmDisplay.leadProbe2.Load sFile, 0, 0
  '      frmDisplay.leadProbe2.DstTop = 0
  '      frmDisplay.leadProbe2.DstLeft = 0
  '      frmDisplay.leadProbe2.DstWidth = frmDisplay.leadProbe2.Width
  '      frmDisplay.leadProbe2.DstHeight = frmDisplay.leadProbe2.Height
  '      frmDisplay.leadProbe2.ForceRepaint
  '      frmDisplay.leadProbe2.Visible = True
    
   ' End If

   ' tmrProbeView.Interval = 100   '100ms
  '  tmrProbeView.Enabled = True

    Exit Sub
    
Error_subDisplayTraining:

    If Err = 20010 Then
        MsgBox Err & ":  " & Error & " - " & sFile
    Else
        MsgBox Err & ":  " & Error
    End If
    
    Exit Sub
    
End Sub

Private Sub subGetFirstTrial()
Dim l As Integer

    
    datSessionResults.Refresh
    datTrialInfo.Recordset.MoveFirst
    For l = 2 To datSessionInfo.Recordset("StartTrial")
        datTrialInfo.Recordset.MoveNext
    Next l
    
    datStimuliInfo.RecordSource = "SELECT * FROM StimuliInfo WHERE TrialNum = " & datSessionInfo.Recordset("StartTrial") & " ORDER BY Position ASC"
    datStimuliInfo.Refresh
    datProbeInfo.RecordSource = "SELECT * FROM ProbeInfo WHERE TrialNum = " & datSessionInfo.Recordset("StartTrial")
    datProbeInfo.Refresh
   
'turn house light on
    With DriverLINXSR1
    .Req_subsystem = DL_DO
    .Req_mode = DL_POLLED
    .Req_op = DL_START
    .Evt_Str_type = DL_NULLEVENT
    .Evt_Stp_type = DL_NULLEVENT
    .Evt_Tim_type = DL_NULLEVENT
    .Sel_chan_format = DL_tNATIVE
    .Sel_chan_N = 1
    .Sel_chan_start = 1   'Channel 0 is configured for Output already
    .Sel_chan_startGainCode = 0
    .Sel_buf_samples = 1
    .Sel_buf_N = 0
    .Refresh
    .Res_Sta_ioValue = APPARATUS_HOUSE_LIGHT
    .Refresh
    Label1.Caption = .Message
End With
    
    tmrHouseLight.Interval = (datTrialInfo.Recordset("ITI") * 1000) - 300
    tmrHouseLight.Enabled = True
    lblStateMessage.Caption = "ITI (" & datTrialInfo.Recordset("ITI") & " sec.)"
    tmrHouseLight1.Interval = datTrialInfo.Recordset("ITI") * 1000
    tmrHouseLight1.Enabled = True
    
    'tmrITI.Interval = datTrialInfo.Recordset("ITI") * 1000
    'flagAbort = False
    'tmrITI.Enabled = True
       
       
    'subInitialItems
    'subTempDisplay
    
    
End Sub

Private Sub subInitTouchscreenRun(tmrInit As Control, sInitState As String, tmrReadTouch As Control)
Stop
 '   tmrInit.Enabled = False
  '  lblStateMessage.Caption = sInitState
   '
'    If sInitState = "Init" Then
  '      agOutp TOUCH_RESET_REGISTER, 0
 '       tmrInit.Interval = 1000
   '     sInitState = "SoftwareReset"
    '    tmrInit.Enabled = True
        
   ' ElseIf sInitState = "SoftwareReset" Then
    '    agOutp TOUCH_DATA_REGISTER, TOUCH_SOFTWARE_RESET
    '    tmrInit.Interval = 1000
    '    sInitState = "ReportTransferOn"
    '    tmrInit.Enabled = True
        
   ' ElseIf sInitState = "ReportTransferOn" Then
   '     agOutp TOUCH_DATA_REGISTER, TOUCH_REPORT_TRANSFER_ON
   '     tmrInit.Interval = 1000
   '     sInitState = "GetErrorReport"
   '     tmrInit.Enabled = True
   '
   ' ElseIf sInitState = "GetErrorReport" Then
   '     agOutp TOUCH_DATA_REGISTER, TOUCH_GET_ERROR_REPORT
   '     tmrInit.Interval = 1000
   '     sInitState = "CoordinateReporting"
   '     tmrInit.Enabled = True
   '
   ' ElseIf sInitState = "CoordinateReporting" Then
   '     agOutp TOUCH_DATA_REGISTER, TOUCH_COORDINATE_REPORTING
   '     tmrInit.Interval = 1000
   '     sInitState = "EnterPointMode"
   '     tmrInit.Enabled = True
   '
   ' ElseIf sInitState = "EnterPointMode" Then
   '     agOutp TOUCH_DATA_REGISTER, TOUCH_ENTER_POINT_MODE
   '     tmrInit.Interval = 1000
   '     sInitState = "StartTrials"
   '     tmrInit.Enabled = True
   '
   ' ElseIf sInitState = "StartTrials" Then
   '     tmrInit.Interval = 0
   '     sInitState = "Init"
   '     tmrInit.Enabled = False
   '
   '     subGetFirstTrial

   ' End If

End Sub

Private Sub subPauseSession()
    
    lblStateMessage.Caption = lblStateMessage.Caption & "[Paused]"
    
    If tmrAbortLTO.Enabled = True Then
        sTimerState = "tmrAbort"
        
    ElseIf tmrIncorrectLTO.Enabled = True Then
        sTimerState = "tmrIncorrectLTO"
        
    ElseIf tmrISI.Enabled = True Then
        sTimerState = "tmrISI"
      
    ElseIf tmrITI.Enabled = True Then
        sTimerState = "tmrITI"
       
    ElseIf tmrProbeDelay.Enabled = True Then
        sTimerState = "tmrProbeDelay"
        
    ElseIf tmrProbeView.Enabled = True Then
        sTimerState = "tmrProbeView"
        
    ElseIf tmrStimuliView.Enabled = True Then
        sTimerState = "tmrStimuliView"
       
    End If
    
    tmrAbortLTO.Enabled = False
    tmrIncorrectLTO.Enabled = False
    tmrISI.Enabled = False
    tmrITI.Enabled = False
    tmrProbeDelay.Enabled = False
    tmrProbeView.Enabled = False
    tmrStimuliView.Enabled = False
    
    mnuStart.Enabled = False
    mnuPause.Enabled = False
    mnuResume.Enabled = True
    mnuReset.Enabled = True
    
End Sub

Private Sub subPercentCorrect()

    Dim iDenominator As Integer
    
    iDenominator = CInt(lblSameTrialSameResp.Caption) _
                    + CInt(lblSameTrialDiffResp.Caption) _
                    + CInt(lblDiffTrialSameResp.Caption) _
                    + CInt(lblDiffTrialDiffResp.Caption)

    
    If iDenominator = 0 Then
        iDenominator = 1
    End If
    
    lblPercentCorrect.Caption = Format(100 * (CInt(lblSameTrialSameResp.Caption) _
                                + CInt(lblDiffTrialDiffResp.Caption)) _
                                / iDenominator, "##0")
    
    datSessionResults.Recordset.MoveFirst
    datSessionResults.Recordset.Edit
    datSessionResults.Recordset("SameCorrect") = lblSameTrialSameResp.Caption
    datSessionResults.Recordset("DiffCorrect") = lblDiffTrialDiffResp.Caption
    datSessionResults.Recordset("SameIncorrect") = lblSameTrialDiffResp.Caption
    datSessionResults.Recordset("DiffIncorrect") = lblDiffTrialSameResp.Caption
    datSessionResults.Recordset.Update

End Sub

Private Sub subResetProbe()
        
    ' mtgsubDisableTouchscreen tmrReadTouch
        
    tmrProbeView.Enabled = False
    tmrProbeView.Interval = 0
    
    'tmrSampleDuration.Enabled = True
    'tmrSampleDuration.Interval = 1  'ms
    
If TrialStatus = True Then
    tmrSampleDuration.Enabled = True
    tmrSampleDuration.Interval = 1  'ms
    tmrCorrectItemDuration.Enabled = True
    tmrCorrectItemDuration.Interval = 8000
End If
If TrialStatus = False Then
    tmrSampleDuration.Enabled = True
    tmrSampleDuration.Interval = 2000  'ms
    tmrCorrectItemDuration.Enabled = True
    tmrCorrectItemDuration.Interval = 2000
End If
 
 
 '   leadStimuli.Visible = False
 '   leadProbe.Visible = False
 '   leadProbe2.Visible = False
   
   ' lblDiffSelection.Visible = False
        
'    subClearDisplay
    
    

End Sub

Private Sub subResetSession()
    
    lblStateMessage.Caption = "Reset"
    
    tmrAbortLTO.Enabled = False
    tmrIncorrectLTO.Enabled = False
    tmrISI.Enabled = False
    tmrITI.Enabled = False
    tmrProbeDelay.Enabled = False
    tmrProbeView.Enabled = False
    tmrStimuliView.Enabled = False
    
    mnuStart.Enabled = True
    mnuPause.Enabled = False
    mnuResume.Enabled = False
    mnuReset.Enabled = False
    
    If datSessionInfo.Recordset("OperationalMode") = "Normal" Then
        mnuTestMode.Enabled = True
        mnuNormalMode.Enabled = False
    Else
        mnuTestMode.Enabled = False
        mnuNormalMode.Enabled = True
    End If
    
    
    leadStimuli.Visible = False
    leadProbe.Visible = False
    leadProbe2.Visible = False
    'lblDiffSelection.Visible = False
    
    subClearDisplay
    
    'clear result tables
    lblSameTrialSameResp.Caption = 0
    lblSameTrialDiffResp.Caption = 0
    lblDiffTrialSameResp.Caption = 0
    lblDiffTrialDiffResp.Caption = 0

End Sub

Private Sub subResumeSession()
    
    If sTimerState = "tmrAbort" Then
        tmrAbortLTO.Enabled = True
        
    ElseIf sTimerState = "tmrIncorrectLTO" Then
        tmrIncorrectLTO.Enabled = True
        
    ElseIf sTimerState = "tmrISI" Then
        tmrISI.Enabled = True
      
    ElseIf sTimerState = "tmrITI" Then
        tmrITI.Enabled = True
       
    ElseIf sTimerState = "tmrProbeDelay" Then
        tmrProbeDelay.Enabled = True
        
    ElseIf sTimerState = "tmrProbeView" Then
        tmrProbeView.Enabled = True
        
    ElseIf sTimerState = "tmrStimuliView" Then
        tmrStimuliView.Enabled = True
       
    End If
    
    sTimerState = ""
    
    mnuStart.Enabled = False
    mnuPause.Enabled = True
    mnuResume.Enabled = False
    mnuReset.Enabled = False
    
End Sub


Private Sub subApparatus()

    Dim sApparatus As String
    
    On Error GoTo Error_subApparatus
    
    waveBeep.Play = 0
    
    sApparatus = datTrialInfo.Recordset("Apparatus")
    lblApparatus.Caption = sApparatus
    
    If sApparatus = "Juice" Then
    
        If CLng(intJuiceDelay.Text) = 0 Then
            tmrApparatusDelay.Interval = 1
        Else
            tmrApparatusDelay.Interval = CLng(intJuiceDelay.Text)
        End If
        
        tmrApparatusDelay.Enabled = True
        
    ElseIf sApparatus = "Pellets" Then
        
        If CLng(intPelletDelay.Text) = 0 Then
            tmrApparatusDelay.Interval = 1
        Else
            tmrApparatusDelay.Interval = CLng(intPelletDelay.Text)
        End If
        
        tmrApparatusDelay.Enabled = True
        
    End If
    
    Exit Sub
    

Error_subApparatus:
    Resume Next
    
End Sub

Private Sub subSessionResults()
    datSessionResults.Refresh
    If datSessionResults.Recordset.RecordCount = 0 Then
        datSessionResults.Recordset.AddNew
    Else
        datSessionResults.Recordset.MoveFirst
        datSessionResults.Recordset.Edit
    End If
    
    datSessionResults.Recordset("SessionID") = lblSessionID.Caption
    datSessionResults.Recordset("Subject") = lblSubject.Caption
    datSessionResults.Recordset("StimTouches") = datSessionInfo.Recordset("StimTouches")
    datSessionResults.Recordset("ProbeTouches") = datSessionInfo.Recordset("ProbeTouches")
    datSessionResults.Recordset("Phase") = datSessionInfo.Recordset("Phase")
    datSessionResults.Recordset("SessionFile") = gsProfile
    datSessionResults.Recordset("Date") = lblDate.Caption
    datSessionResults.Recordset("SetSize") = lblSetSize.Caption
    datSessionResults.Recordset("Sets") = datSessionInfo.Recordset("Sets")
    datSessionResults.Recordset("StartTime") = Now
    datSessionResults.Recordset("EndTime") = Now
    datSessionResults.Recordset("NumOfTrials") = lblNumOfTrials.Caption
    datSessionResults.Recordset("AbortCorrectionEnabled") = chkAbortCorrectionEnabled.Value
    datSessionResults.Recordset("IncorrectCorrectionEnabled") = chkIncorrectCorrectionEnabled.Value
    datSessionResults.Recordset("TotalAborts") = lblNumAborts.Caption
    datSessionResults.Recordset("TotalAbortCorrections") = lblTotalAbortCorrections.Caption
    datSessionResults.Recordset("TotalIncorrectCorrections") = lblTotalIncorrectCorrections.Caption
    datSessionResults.Recordset("PercentCorrect") = 0 'lblPercentCorrect.Caption
    datSessionResults.Recordset("SameCorrect") = 0 'lblSameTrialSameResp.Caption
    datSessionResults.Recordset("DiffCorrect") = 0 'lblDiffTrialDiffResp.Caption
    datSessionResults.Recordset("SameIncorrect") = 0 'lblSameTrialDiffResp.Caption
    datSessionResults.Recordset("DiffIncorrect") = 0 'lblDiffTrialSameResp.Caption
    datSessionResults.Recordset("JuiceProbability") = lblJuiceProbability.Caption
    datSessionResults.Recordset("PelletProbability") = lblPelletProbability.Caption
    datSessionResults.Recordset("ReinfProbability") = lblReinfProbability.Caption
    datSessionResults.Recordset("SameProbability") = lblSameProbability.Caption
    datSessionResults.Recordset("DiffProbability") = lblDiffProbability.Caption
   ' datSessionResults.Recordset("VariableListLength") = chkVariableListLength.Value
    datSessionResults.Recordset.Update

End Sub
Private Sub subSetTestMode()
        
    mnuNormalMode.Enabled = True
    mnuTestMode.Enabled = False

    lblMode.BackColor = COLOR_RED
    lblMode.ForeColor = COLOR_WHITE
    
    If datSessionInfo.Recordset("TrainingMode") = "Off" Then
        lblMode.Caption = "Test Mode"
    Else
        lblMode.Caption = "Training - Test Mode"
    End If
    
    datSessionInfo.Recordset.Edit
    datSessionInfo.Recordset("OperationalMode") = "Test"
    datSessionInfo.Recordset.Update
    
'mt    gsubDisableTouchscreen tmrReadTouch
   ' agOutp METRABYTE_CONTROL, APPARATUS_SET_B_OUTPUT 'APPARATUS_SET_B_INPUT

End Sub

Private Sub subAbort()

    lblNumAborts.Caption = lblNumAborts.Caption + 1
    
    If datSessionInfo.Recordset("ProbeAbortLTO") > 0 Then
        lblStateMessage.Caption = "ProbeAbortLTO (" & datSessionInfo.Recordset("ProbeAbortLTO") & " msec.)"
        tmrProbeAbort.Interval = datSessionInfo.Recordset("ProbeAbortLTO")
        tmrProbeAbort.Enabled = True
    Else
        lblStateMessage.Caption = "ITI (" & datTrialInfo.Recordset("ITI") & " sec.)"
        tmrITI.Interval = datTrialInfo.Recordset("ITI") * 1000
        tmrITI.Enabled = True
    End If
End Sub

Private Sub subAbortLTO()

    'modified for aborts through out trial
    tmrISI.Enabled = False
    tmrProbeDelay.Enabled = False
    tmrProbeView.Enabled = False
    tmrStimuliView.Enabled = False
    flagAbort = False
    Do
       ' DoEvents
        flagAbort = True
    Loop Until tmrITI.Enabled = False And tmrIncorrectLTO.Enabled = False
       
    lblStateMessage.Caption = "Abort LTO/Correction (" & datSessionInfo.Recordset("DefaultAbortLTO") & " sec.)"
    lblNumAborts.Caption = lblNumAborts.Caption + 1
    lblTotalAbortCorrections.Caption = lblTotalAbortCorrections.Caption + 1
    tmrAbortLTO.Interval = datSessionInfo.Recordset("DefaultAbortLTO") * 1000
    tmrAbortLTO.Enabled = True
    
End Sub


Private Sub subCorrect()

    subApparatus
    If datTrialInfo.Recordset("Transfer") = False Then
        subPercentCorrect
    Else
        subTranPercentCorrect
    End If
      
End Sub
    
    
Private Sub subInitialItems()
    Dim ttype As Integer
    ttype = datTrialInfo.Recordset("Configuration")

    Randomize
    Dim lViewTime As Long
    Dim trip As Integer
    
    On Error GoTo Error_subInitialItems
    
    If datTrialInfo.Recordset("NoTouchTrial") = -1 Then  'ts fix ?
         sResponse = "None"
    End If
        
    RespLoc = 0
    lResponseTime = 0
    sImageDisplayed = "first"
    gPecks = 0
    'pos1
    stimulusFiles(1) = datSessionInfo.Recordset("ImagePath") & "\" & datStimuliInfo.Recordset("Filename")
    'pos2
    stimulusFiles(2) = datSessionInfo.Recordset("ImagePath") & "\" & datProbeInfo.Recordset("Filename")
    
    'for the rest of the stims we can easily loop over them
    Dim stim As Integer
    For stim = 3 To ttype
        stimulusFiles(stim) = datSessionInfo.Recordset("ImagePath") & "\" _
            & datProbeInfo.Recordset("Filename" & (stim - 1))
    Next stim
    
    stimulusFiles(PROBE_LOCATION) = datSessionInfo.Recordset("ImagePath") _
                    & "\" & datProbeInfo.Recordset("ProbeFileName")
    
    lViewTime = datStimuliInfo.Recordset("Viewtime")
    gDisplayTime = lViewTime
    lblStateMessage.Caption = "Initial (" & lViewTime & " sec.): "
    
    With leadStimuli
        .AutoScroll = False
        .Load stimulusFiles(1), 0, 0
        .DstTop = 0
        .DstLeft = 0
        .DstWidth = leadStimuli.Width / 2
        .DstHeight = leadStimuli.Height / 2
        .ForceRepaint
        .Visible = True
    End With
    With leadProbe
        .AutoScroll = False
        .Load stimulusFiles(2), 0, 0
        .DstTop = 0
        .DstLeft = 0
        .DstWidth = leadProbe.Width / 2
        .DstHeight = leadProbe.Height / 2
        .ForceRepaint
        .Visible = True
    End With
    
    Dim p As Integer
    For p = 1 To ttype
        positions(p) = Int(datTrialInfo.Recordset("POS" & p))
    Next p
    
    Dim currentStim As Integer
    
    'loop throught the appropriate number of stimuli and display the image
    For currentStim = 1 To ttype
        'note that we are making use of the positions and stimulusFiles array
        With stimulus_array(positions(currentStim))
            .Load stimulusFiles(currentStim), 0, 0
            .DstTop = 0
            .DstLeft = 0
            .DstWidth = .Width
            .DstHeight = .Height
            '.Flip
            '.Reverse
            .ForceRepaint
            .Visible = True
        End With
    Next
    
    tmrProbeView.Interval = 100   '100ms
    tmrProbeView.Enabled = True

    Exit Sub
    
Error_subInitialItems:

    If Err = 20010 Then
        MsgBox Err & ":  " & Error & " - " & sFile
    Else
        MsgBox Err & ":  " & Error
    End If
    
    Exit Sub
    
End Sub
Private Sub subSecondItems()

    Dim lViewTime As Long
   
    On Error GoTo Error_subSecondItems
    
    If datTrialInfo.Recordset("NoTouchTrial") = -1 Then  'ts fix ?
         sResponse = "None"
    End If
    
    RespLoc = 0
    lResponseTime = 0
    sImageDisplayed = "second"
    lblStateMessage.Caption = "Second (" & lViewTime & " sec.): "
     
    leadStimuli.AutoScroll = False
    leadStimuli.Load stimulusFiles(1), 0, 0
    leadStimuli.DstTop = 0
    leadStimuli.DstLeft = 0
    leadStimuli.DstWidth = leadStimuli.Width / 2
    leadStimuli.DstHeight = leadStimuli.Height / 2
    leadStimuli.ForceRepaint
    leadStimuli.Visible = True
    
    leadProbe.AutoScroll = False
    leadProbe.Load stimulusFiles(2), 0, 0
    leadProbe.DstTop = 0
    leadProbe.DstLeft = 0
    leadProbe.DstWidth = leadProbe.Width / 2
    leadProbe.DstHeight = leadProbe.Height / 2
    leadProbe.ForceRepaint
    leadProbe.Visible = True
    
    leadProbe2.AutoScroll = False
    leadProbe2.Load stimulusFiles(PROBE_LOCATION), 0, 0
    leadProbe2.DstTop = 0
    leadProbe2.DstLeft = 0
    leadProbe2.DstWidth = leadProbe2.Width / 2
    leadProbe2.DstHeight = leadProbe2.Height / 2
    leadProbe2.ForceRepaint
    leadProbe2.Visible = True
    
If datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1") Then
   With stimulus_array(positions(1))
        .Load stimulusFiles(PROBE_LOCATION), 0, 0 'as mentioned at the top of the file, 7 is always probe
        .DstTop = 0
        .DstLeft = 0
        .DstWidth = .Width
        .DstHeight = .Height
        '.Flip
        '.Reverse
        .ForceRepaint
        .Visible = True
   End With
    
   With stimulus_array(positions(2))
        .Load stimulusFiles(2), 0, 0  'as mentioned at the top of the file, 2 is always "no change"
        .DstTop = 0
        .DstLeft = 0
        .DstWidth = .Width
        .DstHeight = .Height
        '.Flip
        '.Reverse
        .ForceRepaint
        .Visible = True
    End With
Else
    MsgBox ("Trying to change position other than POS1--Something bad has happened...")
    Stop
End If
    
    
    tmrProbeView.Interval = 100   '100ms
    tmrProbeView.Enabled = True

    Exit Sub
    
Error_subSecondItems:

    If Err = 20010 Then
        MsgBox Err & ":  " & Error & " - " & sFile
    Else
        MsgBox Err & ":  " & Error
    End If
    
    Exit Sub
    
End Sub

Private Sub subTempDisplay()

    Dim sFile As String
    Dim sFile2 As String
    Dim lViewTime As Long
    
    On Error GoTo Error_subTempDisplay
    
    If datTrialInfo.Recordset("NoTouchTrial") = -1 Then  'ts fix ?
         sResponse = "None"
    End If
    
    lResponseTime = 0
    sImageDisplayed = "probe"
    sFile = datSessionInfo.Recordset("ImagePath") & "\" & datProbeInfo.Recordset("Filename")
    ItemP = datProbeInfo.Recordset("Filename")
    sFile2 = datSessionInfo.Recordset("ImagePath") & "\" & datProbeInfo.Recordset("Filename2")
    ItemP2 = datProbeInfo.Recordset("Filename2")
    lblStateMessage.Caption = "Choice (" & lViewTime & " sec.): " & sFile

    
       
    leadStimuli.AutoScroll = False
    leadStimuli.Load sFile, 0, 0
    leadStimuli.DstTop = 0
    leadStimuli.DstLeft = 0
    leadStimuli.DstWidth = leadStimuli.Width
    leadStimuli.DstHeight = leadStimuli.Height
    leadStimuli.ForceRepaint
    leadStimuli.Visible = True
    
    'Lead1.AutoScroll = False
    'Lead1.Load sFile2, 0, 0
    'Lead1.DstTop = 0
    'Lead1.DstLeft = 0
    'Lead1.DstWidth = Lead1.Width
    'Lead1.DstHeight = Lead1.Height
    'Lead1.ForceRepaint
    'Lead1.Visible = True
    
    If datSessionInfo.Recordset("OperationalMode") = "Test" Then
       ' lblDiffSelection.Visible = True
        
    ElseIf datSessionInfo.Recordset("OperationalMode") = "Normal" Then
        frmDisplay.leadStimuli.Load sFile, 0, 0
        'frmDisplay.leadStimuli.DstTop = 0
        'frmDisplay.leadStimuli.DstLeft = 0
        'frmDisplay.leadStimuli.DstWidth = frmDisplay.leadStimuli.Width
        'frmDisplay.leadStimuli.DstHeight = frmDisplay.leadStimuli.Height
        frmDisplay.leadStimuli.Flip
        frmDisplay.leadStimuli.Reverse
        frmDisplay.leadStimuli.ForceRepaint
        frmDisplay.leadStimuli.Visible = True
'        frmDisplay.lblDiff.Visible = True      'takeout shape

        frmDisplay.Lead1.Load sFile2, 0, 0
        'frmDisplay.Lead1.DstTop = 0
        'frmDisplay.Lead1.DstLeft = 0
        'frmDisplay.Lead1.DstWidth = frmDisplay.Lead1.Width
        'frmDisplay.Lead1.DstHeight = frmDisplay.Lead1.Height
        frmDisplay.Lead1.Flip
        frmDisplay.Lead1.Reverse
        frmDisplay.Lead1.ForceRepaint
        frmDisplay.Lead1.Visible = True

        frmDisplay.Lead2.Load sFile2, 0, 0
        frmDisplay.Lead2.Flip
        frmDisplay.Lead2.Reverse
        frmDisplay.Lead2.ForceRepaint
        frmDisplay.Lead2.Visible = True
        
        frmDisplay.Lead3.Load sFile2, 0, 0
        frmDisplay.Lead3.Flip
        frmDisplay.Lead3.Reverse
        frmDisplay.Lead3.ForceRepaint
        frmDisplay.Lead3.Visible = True

        frmDisplay.Lead4.Load sFile2, 0, 0
        frmDisplay.Lead4.Flip
        frmDisplay.Lead4.Reverse
        frmDisplay.Lead4.ForceRepaint
        frmDisplay.Lead4.Visible = True

        frmDisplay.Lead5.Load sFile2, 0, 0
        frmDisplay.Lead5.Flip
        frmDisplay.Lead5.Reverse
        frmDisplay.Lead5.ForceRepaint
        frmDisplay.Lead5.Visible = True

        frmDisplay.Lead6.Load sFile2, 0, 0
        frmDisplay.Lead6.Flip
        frmDisplay.Lead6.Reverse
        frmDisplay.Lead6.ForceRepaint
        frmDisplay.Lead6.Visible = True

        frmDisplay.Lead7.Load sFile2, 0, 0
        frmDisplay.Lead7.Flip
        frmDisplay.Lead7.Reverse
        frmDisplay.Lead7.ForceRepaint
        frmDisplay.Lead7.Visible = True

        frmDisplay.Lead8.Load sFile2, 0, 0
        frmDisplay.Lead8.Flip
        frmDisplay.Lead8.Reverse
        frmDisplay.Lead8.ForceRepaint
        frmDisplay.Lead8.Visible = True

        frmDisplay.Lead9.Load sFile2, 0, 0
        frmDisplay.Lead9.Flip
        frmDisplay.Lead9.Reverse
        frmDisplay.Lead9.ForceRepaint
        frmDisplay.Lead9.Visible = True

        frmDisplay.Lead10.Load sFile2, 0, 0
        frmDisplay.Lead10.Flip
        frmDisplay.Lead10.Reverse
        frmDisplay.Lead10.ForceRepaint
        frmDisplay.Lead10.Visible = True

        frmDisplay.Lead11.Load sFile2, 0, 0
        frmDisplay.Lead11.Flip
        frmDisplay.Lead11.Reverse
        frmDisplay.Lead11.ForceRepaint
        frmDisplay.Lead11.Visible = True
        
        frmDisplay.Lead12.Load sFile2, 0, 0
        frmDisplay.Lead12.Flip
        frmDisplay.Lead12.Reverse
        frmDisplay.Lead12.ForceRepaint
        frmDisplay.Lead12.Visible = True

        frmDisplay.Lead13.Load sFile2, 0, 0
        frmDisplay.Lead13.Flip
        frmDisplay.Lead13.Reverse
        frmDisplay.Lead13.ForceRepaint
        frmDisplay.Lead13.Visible = True

        frmDisplay.Lead14.Load sFile2, 0, 0
        frmDisplay.Lead14.Flip
        frmDisplay.Lead14.Reverse
        frmDisplay.Lead14.ForceRepaint
        frmDisplay.Lead14.Visible = True

        frmDisplay.Lead15.Load sFile2, 0, 0
        frmDisplay.Lead15.Flip
        frmDisplay.Lead15.Reverse
        frmDisplay.Lead15.ForceRepaint
        frmDisplay.Lead15.Visible = True

'        If datTrialInfo.Recordset("TrialType") = "Same" Then frmDisplay.lblProbe.Visible = True  'takeout shape
'        If datTrialInfo.Recordset("TrialType") = "Diff" Then frmDisplay.lblDiff.Visible = True  'takeout shape
     
    End If
    
    tmrProbeView.Interval = 100   '100ms
    tmrProbeView.Enabled = True

    Exit Sub
    
Error_subTempDisplay:

    If Err = 20010 Then
        MsgBox Err & ":  " & Error & " - " & sFile
    Else
        MsgBox Err & ":  " & Error
    End If
    
    Exit Sub
    

End Sub



Private Sub subDisplayProbe()

    Dim sFile As String
    Dim sFile2 As String
    Dim lViewTime As Long
    
   
   ' On Error GoTo Error_subDisplayProbe
    
 '   If datTrialInfo.Recordset("NoTouchTrial") = -1 Then  'ts fix ?
 '        sResponse = "None"
 '   End If
    
 '   lResponseTime = 0
 '   sImageDisplayed = "probe"
 '   sFile = datSessionInfo.Recordset("ImagePath") & "\" & datProbeInfo.Recordset("Filename")
 '   ItemP = datProbeInfo.Recordset("Filename")
 '   sFile2 = datSessionInfo.Recordset("ImagePath") & "\" & datProbeInfo.Recordset("Filename2")
 '   ItemP2 = datProbeInfo.Recordset("Filename2")
 '   lblStateMessage.Caption = "Choice (" & lViewTime & " sec.): " & sFile

 '   leadProbe.AutoScroll = False
 '   leadProbe.Load sFile, 0, 0
 '   leadProbe.DstTop = 0
 '   leadProbe.DstLeft = 0
 '   leadProbe.DstWidth = leadProbe.Width
 '   leadProbe.DstHeight = leadProbe.Height
 '   leadProbe.ForceRepaint
 '   leadProbe.Visible = True
 '
 '   leadProbe2.AutoScroll = False
 '   leadProbe2.Load sFile2, 0, 0
 '   leadProbe2.DstTop = 0
 '   leadProbe2.DstLeft = 0
 '   leadProbe2.DstWidth = leadProbe2.Width
 '   leadProbe2.DstHeight = leadProbe2.Height
 '   leadProbe2.ForceRepaint
 '   leadProbe2.Visible = True
 '
 'if datSessionInfo.Recordset("OperationalMode") = "Test" Then
       ' lblDiffSelection.Visible = True
        
  '  ElseIf datSessionInfo.Recordset("OperationalMode") = "Normal" Then
  '      frmDisplay.leadProbe.Load sFile, 0, 0
  ' '     frmDisplay.leadProbe.DstTop = 0
  ' '     frmDisplay.leadProbe.DstLeft = 0
  ' '     frmDisplay.leadProbe.DstWidth = frmDisplay.leadProbe.Width
  ' '     frmDisplay.leadProbe.DstHeight = frmDisplay.leadProbe.Height
  '      frmDisplay.leadProbe.Flip
  '      frmDisplay.leadProbe.Reverse
  '      frmDisplay.leadProbe.ForceRepaint
  '      frmDisplay.leadProbe.Visible = True    'takeout shape
' '       frmDisplay.lblDiff.Visible = True      'takeout shape'

        'frmDisplay.leadProbe2.Load sFile2, 0, 0
    '    frmDisplay.leadProbe2.DstTop = 0
    '    frmDisplay.leadProbe2.DstLeft = 0
    '    frmDisplay.leadProbe2.DstWidth = frmDisplay.leadProbe2.Width
    '    frmDisplay.leadProbe2.DstHeight = frmDisplay.leadProbe2.Height
    '    frmDisplay.leadProbe2.Flip
    '    frmDisplay.leadProbe2.Reverse
    '    frmDisplay.leadProbe2.ForceRepaint
    '    frmDisplay.leadProbe2.Visible = True    'takeout shape

'        If datTrialInfo.Recordset("TrialType") = "Same" Then frmDisplay.lblProbe.Visible = True  'takeout shape
'        If datTrialInfo.Recordset("TrialType") = "Diff" Then frmDisplay.lblDiff.Visible = True  'takeout shape
     
   ' End If
    
   ' tmrProbeView.Interval = 100   '100ms
   ' tmrProbeView.Enabled = True

    Exit Sub
    
Error_subDisplayProbe:

    If Err = 20010 Then
        MsgBox Err & ":  " & Error & " - " & sFile
    Else
        MsgBox Err & ":  " & Error
    End If
    
    Exit Sub
    

End Sub

Private Sub subDisplayStimuli()

    Dim sFile As String
    Dim lViewTime As Long
    
    On Error GoTo Error_subDisplayStimuli
    cursor_visible (False)
      
    sImageDisplayed = "stim"
    iStimTouches = 0
    lResponseTime = 0
    sResponse = "None"
    lblResponse.Caption = ""
    lblApparatus.Caption = ""
    
    lblISI.Caption = datStimuliInfo.Recordset("ISI")
    
    lViewTime = datStimuliInfo.Recordset("ViewTime")
    
    sFile = datSessionInfo.Recordset("ImagePath") & "\" & datStimuliInfo.Recordset("Filename")
    
    subRecordItem
    
    lblStateMessage.Caption = "Stimuli " & datStimuliInfo.Recordset("Position") & " (" & lViewTime & " msec.): " & sFile
    lblViewTime.Caption = lViewTime
    
    leadStimuli.AutoScroll = False
    leadStimuli.Load sFile, 0, 0
    leadStimuli.DstTop = 0
    leadStimuli.DstLeft = 0
    leadStimuli.DstWidth = leadStimuli.Width
    leadStimuli.DstHeight = leadStimuli.Height
    leadStimuli.ForceRepaint
    leadStimuli.Visible = True  'takeout shape
    leadProbe.Visible = False
    leadProbe2.Visible = False
  '  lblDiffSelection.Visible = False
    
    If datSessionInfo.Recordset("OperationalMode") = "Normal" Then
        frmDisplay.leadStimuli.Load sFile, 0, 0
        'frmDisplay.leadStimuli.DstTop = 10
        'frmDisplay.leadStimuli.DstLeft = 0
        'frmDisplay.leadStimuli.DstWidth = frmDisplay.leadStimuli.Width
        'frmDisplay.leadStimuli.DstHeight = frmDisplay.leadStimuli.Height
        frmDisplay.leadStimuli.Flip
        frmDisplay.leadStimuli.Reverse
        frmDisplay.leadStimuli.ForceRepaint
        frmDisplay.leadStimuli.Visible = True 'takeout shape
    End If
    
    If datTrialInfo.Recordset("NoTouchTrial") = -1 Then
        'True = no touch trial
        tmrStimuliView.Interval = lViewTime '* 1000
        tmrStimuliView.Enabled = True
    End If
    
    tmrProbeView.Interval = 100   '100ms
    tmrProbeView.Enabled = True
    
    Exit Sub
    
Error_subDisplayStimuli:

    If Err = 20010 Then
        MsgBox Err & ":  " & Error & " - " & sFile
    Else
        MsgBox Err & ":  " & Error
    End If
    
    Exit Sub
    
End Sub


Private Sub subIncorrect()

    If datTrialInfo.Recordset("Transfer") = False Then
        subPercentCorrect
    Else
        subTranPercentCorrect
    End If
    
   'turn house light on
    With DriverLINXSR1
    .Req_subsystem = DL_DO
    .Req_mode = DL_POLLED
    .Req_op = DL_START
    .Evt_Str_type = DL_NULLEVENT
    .Evt_Stp_type = DL_NULLEVENT
    .Evt_Tim_type = DL_NULLEVENT
    .Sel_chan_format = DL_tNATIVE
    .Sel_chan_N = 1
    .Sel_chan_start = 1   'Channel 0 is configured for Output already
    .Sel_chan_startGainCode = 0
    .Sel_buf_samples = 1
    .Sel_buf_N = 0
    .Refresh
    .Res_Sta_ioValue = APPARATUS_HOUSE_LIGHT
    .Refresh
    Label1.Caption = .Message
End With
    
    tmrHouseLight.Interval = (datTrialInfo.Recordset("ITI") * 1000) - 300
    tmrHouseLight.Enabled = True
    
   ' If datSessionInfo.Recordset("IncorrectCorrectionsEnabled") Then
       ' tmrApparatusDelay.Interval = 1
       ' tmrApparatusDelay.Enabled = True
      '  lblStateMessage.Caption = "Time Out (" & datSessionInfo.Recordset("DefaultIncorrectLTO") & " sec.)"
      '  lblTotalIncorrectCorrections.Caption = lblTotalIncorrectCorrections.Caption + 1
      '  tmrIncorrectLTO.Interval = datSessionInfo.Recordset("DefaultIncorrectLTO") * 1000
      '  tmrIncorrectLTO.Enabled = True
'    Else
    lblStateMessage.Caption = "ITI (" & datTrialInfo.Recordset("ITI") & " sec.)"
    tmrITI.Interval = datTrialInfo.Recordset("ITI") * 1000
    'End If
   flagAbort = False
   tmrITI.Enabled = True
    
End Sub

Private Sub subIncorrectLTO()
    
    If datTrialInfo.Recordset("Transfer") = False Then
        subPercentCorrect
    Else
        subTranPercentCorrect
    End If
    lblStateMessage.Caption = "Time Out (" & datSessionInfo.Recordset("DefaultIncorrectLTO") & " sec.)"
    lblTotalIncorrectCorrections.Caption = lblTotalIncorrectCorrections.Caption + 1
    
    tmrIncorrectLTO.Interval = datSessionInfo.Recordset("DefaultIncorrectLTO") * 1000
    tmrIncorrectLTO.Enabled = True
           
End Sub

Private Sub subSetNormalMode()
    
    mnuNormalMode.Enabled = False
    mnuTestMode.Enabled = True

    lblMode.BackColor = COLOR_GREEN
    lblMode.ForeColor = COLOR_BLACK
    
    If datSessionInfo.Recordset("TrainingMode") = "Off" Then
        lblMode.Caption = "Normal Mode"
    Else
        lblMode.Caption = "Training - Normal Mode"
    End If
    
    datSessionInfo.Recordset.Edit
    datSessionInfo.Recordset("OperationalMode") = "Normal"
    datSessionInfo.Recordset.Update
    
'mt    gsubEnableTouchscreen tmrReadTouch
    'agOutp METRABYTE_CONTROL, APPARATUS_SET_B_OUTPUT
    'initialize PCI driver
    With DriverLINXSR1
        .Req_device = 0
        .Req_mode = DL_OTHER
        .Req_op = DL_INITIALIZE
        .Req_DLL_name = "KPCIPIO"
        .Refresh
        Label1.Caption = .Message
    End With
    
    ' code below will Configure Channel 1 for Output
    With DriverLINXSR1
        .Req_subsystem = DL_DO  ' the subsystem to assign the channel to
        .Req_mode = DL_OTHER
        .Req_op = DL_CONFIGURE  ' it is a configuration type operation
        .Evt_Tim_type = DL_DIOSETUP
        .Evt_Tim_dioChannel = 1   ' configure channel 1 or Port b as Output
        .Evt_Tim_dioMode = DL_DIO_BASIC
        .Evt_Str_type = DL_NULLEVENT
        .Evt_Stp_type = DL_NULLEVENT
        .Sel_chan_N = 0
        .Refresh
        lblStateMessage.Caption = .Message
    End With
    
    
End Sub


Private Sub subStartSession()
    
    lblStateMessage.Caption = "Start"
    
    mnuNormalMode.Enabled = False
    mnuTestMode.Enabled = False
    
    mnuStart.Enabled = False
    mnuPause.Enabled = True
    mnuResume.Enabled = False
    mnuReset.Enabled = False
    
 
    'reset response grid
    If datSessionInfo.Recordset("StartTrial") = 1 Then
        lblSameTrialSameResp.Caption = 0
        lblSameTrialDiffResp.Caption = 0
        lblDiffTrialSameResp.Caption = 0
        lblDiffTrialDiffResp.Caption = 0
    
        lblPercentCorrect.Caption = "0"
        lblNumAborts.Caption = "0"
        lblNumJuice.Caption = "0"
        lblNumPellets.Caption = "0"
        lblTotalAbortCorrections.Caption = "0"
        lblTotalIncorrectCorrections.Caption = "0"
    End If
    
    lResponseTime = 0
    sResponse = "None"
    lblResponse.Caption = ""
    lblApparatus.Caption = ""
    

     If datSessionInfo.Recordset("StartTrial") = 1 Then subSessionResults
    
    If datSessionInfo.Recordset("OperationalMode") = "Normal" Then
        'agOutp METRABYTE_CONTROL, APPARATUS_SET_B_OUTPUT
        
        frmDisplay.Show
        cursor_visible (False)
        subGetFirstTrial
    
    ElseIf datSessionInfo.Recordset("OperationalMode") = "Test" Then

        'agOutp METRABYTE_CONTROL, APPARATUS_SET_B_OUTPUT
        subGetFirstTrial
    
    End If
    
End Sub

Private Sub subStopSession()

    Dim iTotalAbortCorrections As Integer
    Dim iTotalIncorrectCorrections As Integer
    
    lblStateMessage.Caption = "Complete"
    
    tmrAbortLTO.Enabled = False
    tmrIncorrectLTO.Enabled = False
    tmrISI.Enabled = False
    tmrITI.Enabled = False
    tmrProbeDelay.Enabled = False
    tmrProbeView.Enabled = False
    tmrStimuliView.Enabled = False
    
    Unload frmDisplay
    
    'Session Summary
    datRecordCounter.DatabaseName = gsResultsFile
    datRecordCounter.RecordSource = "SELECT SUM(AbortCorrections)AS Quantity FROM TrialResults"
    datRecordCounter.Refresh
    iTotalAbortCorrections = datRecordCounter.Recordset("Quantity")
    
    datRecordCounter.DatabaseName = gsResultsFile
    datRecordCounter.RecordSource = "SELECT SUM(IncorrectCorrections) as Quantity FROM TrialResults"
    datRecordCounter.Refresh
    iTotalIncorrectCorrections = datRecordCounter.Recordset("Quantity")
    
    datSessionResults.Recordset.MoveFirst
    datSessionResults.Recordset.Edit
    datSessionResults.Recordset("EndTime") = Now
    datSessionResults.Recordset("AbortCorrectionEnabled") = chkAbortCorrectionEnabled.Value
    datSessionResults.Recordset("IncorrectCorrectionEnabled") = chkIncorrectCorrectionEnabled.Value
    datSessionResults.Recordset("TotalAborts") = lblNumAborts.Caption
    datSessionResults.Recordset("TotalAbortCorrections") = iTotalAbortCorrections
    datSessionResults.Recordset("TotalIncorrectCorrections") = iTotalIncorrectCorrections
    datSessionResults.Recordset("PercentCorrect") = lblPercentCorrect.Caption
    datSessionResults.Recordset("SameCorrect") = lblSameTrialSameResp.Caption
    datSessionResults.Recordset("DiffCorrect") = lblDiffTrialDiffResp.Caption
    datSessionResults.Recordset("SameIncorrect") = lblSameTrialDiffResp.Caption
    datSessionResults.Recordset("DiffIncorrect") = lblDiffTrialSameResp.Caption
    datSessionResults.Recordset.Update
    
    mnuStart.Enabled = False
    mnuPause.Enabled = False
    mnuResume.Enabled = False
    mnuReset.Enabled = True
    
    leadStimuli.Visible = False
    leadProbe.Visible = False
    leadProbe2.Visible = False
'    lblDiffSelection.Visible = False
    
    cursor_visible (True)
    
End Sub

Private Sub subTrialResults()
        
    datTrialResults.Recordset.AddNew
    datTrialResults.Recordset("TrialNum") = lblTrialNum.Caption
    datTrialResults.Recordset("Configuration") = datTrialInfo.Recordset("Configuration")
    datTrialResults.Recordset("TrialType") = lblType.Caption
    datTrialResults.Recordset("Apparatus") = datTrialInfo.Recordset("Apparatus")
    datTrialResults.Recordset("Response") = RespLoc 'HoldResponse 'sResponse
    datTrialResults.Recordset("ResponseTime") = lResponseTime
    datTrialResults.Recordset("CumPercentCorrect") = lblPercentCorrect.Caption
    datTrialResults.Recordset("ListLength") = lblListLength.Caption
    datTrialResults.Recordset("ProbeDelay") = lblProbeDelay.Caption
    
    datTrialResults.Recordset("ITI") = lblITI.Caption
    
    datTrialResults.Recordset("CorrectPosition") = datTrialInfo.Recordset("POS1")
    
    datTrialResults.Recordset("AbortCorrections") = lblTotalAbortCorrections.Caption
    datTrialResults.Recordset("IncorrectCorrections") = lblTotalIncorrectCorrections.Caption
    
    datTrialResults.Recordset("Item1") = datStimuliInfo.Recordset("Filename")
    datTrialResults.Recordset("Item2") = datProbeInfo.Recordset("Filename")
    datTrialResults.Recordset("Item3") = datProbeInfo.Recordset("Filename2")
    datTrialResults.Recordset("Item4") = datProbeInfo.Recordset("Filename3")
    datTrialResults.Recordset("Item5") = datProbeInfo.Recordset("Filename4")
    datTrialResults.Recordset("Item6") = datProbeInfo.Recordset("Filename5")
    datTrialResults.Recordset("Item7") = datProbeInfo.Recordset("Filename6")
    datTrialResults.Recordset("Item8") = datProbeInfo.Recordset("Filename7")
    
    datTrialResults.Recordset("ItemP") = datProbeInfo.Recordset("Filename")
    datTrialResults.Recordset("ItemP2") = datProbeInfo.Recordset("ProbeFileName")
    
    datTrialResults.Recordset("Pos1") = datTrialInfo.Recordset("POS1")
    datTrialResults.Recordset("Pos2") = datTrialInfo.Recordset("POS2")
    datTrialResults.Recordset("Pos3") = datTrialInfo.Recordset("POS3")
    datTrialResults.Recordset("Pos4") = datTrialInfo.Recordset("POS4")
    datTrialResults.Recordset("Pos5") = datTrialInfo.Recordset("POS5")
    datTrialResults.Recordset("Pos6") = datTrialInfo.Recordset("POS6")
    datTrialResults.Recordset("Pos7") = datTrialInfo.Recordset("POS7")
    datTrialResults.Recordset("Pos8") = datTrialInfo.Recordset("POS8")
    
    datTrialResults.Recordset("Transfer") = datTrialInfo.Recordset("Transfer")
    datTrialResults.Recordset("NoTouch") = datTrialInfo.Recordset("NoTouchTrial")
    datTrialResults.Recordset.Update
        
    lblTotalAbortCorrections.Caption = 0
    lblTotalIncorrectCorrections.Caption = 0

End Sub

Private Sub Form_Load()

    On Error GoTo Error_Form_Load

    leadStimuli.BackColor = 0
    leadProbe.BackColor = 0
    leadProbe2.BackColor = 0
    
    gsubScreenDimensions Me, 0, 0
    
    mnuPause.Enabled = False
    mnuResume.Enabled = False
    lblDate.Caption = Format(Now, "yy/mm/dd hh:mm:ss")
    
    lblApparatus.Caption = ""
    lblResponse.Caption = ""
    
    Refresh
    datSessionResults.DatabaseName = gsResultsFile
    datSessionResults.RecordSource = "SessionResults"
    datSessionResults.Refresh
    
    datTrialResults.DatabaseName = gsResultsFile
    datTrialResults.RecordSource = "TrialResults"
    datTrialResults.Refresh
    
    datStimResponses.DatabaseName = gsResultsFile
    datStimResponses.RecordSource = "StimResponses"
    datStimResponses.Refresh
    
    datSessionInfo.DatabaseName = gsProfile
    datSessionInfo.RecordSource = "SessionInfo"
    datSessionInfo.Refresh

    waveBeep.WaveIndex = 0
    waveBeep.FilenameRead = datSessionInfo.Recordset("SoundFile")
    
    datTrialInfo.DatabaseName = gsProfile
    datTrialInfo.RecordSource = "TrialInfo"
    datTrialInfo.Refresh

    datStimuliInfo.DatabaseName = gsProfile
    datProbeInfo.DatabaseName = gsProfile

    If datSessionInfo.Recordset("OperationalMode") = "Normal" Then
        subSetNormalMode
    ElseIf datSessionInfo.Recordset("OperationalMode") = "Test" Then
        subSetTestMode
    End If
    
    Exit Sub
    

Error_Form_Load:

    Resume Next
    
End Sub



Private Sub Form_MouseDown(Button As Integer, Shift As Integer, X As Single, Y As Single)
butt = Button
End Sub

Private Sub lblDiffSelection_Click()

    If datSessionInfo.Recordset("OperationalMode") = "Test" Then
    
        sResponse = "Diff"
        subCheckResponse

    End If

End Sub

Private Sub leadProbe_Click()

    If datSessionInfo.Recordset("OperationalMode") = "Test" Then
    
        sResponse = "Same"
        subCheckResponse

    End If
    
End Sub

Private Sub leadStimuli_Click()

    If datSessionInfo.Recordset("OperationalMode") = "Test" Then
    
        sResponse = "Stim"
        subCheckResponse

    End If

End Sub

Private Sub mnuClose_Click()
    
    If datSessionInfo.Recordset("OperationalMode") = "Normal" Then
        Unload frmDisplay
    End If
    
    cursor_visible (True)
    Unload Me
    
End Sub


Private Sub mnuNormalMode_Click()

    subSetNormalMode
    
End Sub

Private Sub mnuPause_Click()

    subPauseSession
    
End Sub

Private Sub mnuReset_Click()

    subResetSession
    
End Sub


Private Sub mnuResume_Click()

    subResumeSession
    
End Sub


Private Sub mnuStart_Click()

    subStartSession

End Sub

Private Sub mnuTestMode_Click()
    
    subSetTestMode
    
End Sub



Private Sub tmrAbortLTO_Timer()

    tmrAbortLTO.Enabled = False
    tmrAbortLTO.Interval = 0
    flagAbort = True
    datStimuliInfo.Recordset.MoveFirst

    'If datSessionInfo.Recordset("TrainingMode") = "Off" Then
     '   subDisplayStimuli
    'Else
    '    subDisplayTraining
    'End If
    
    lblStateMessage.Caption = "ITI (" & datTrialInfo.Recordset("ITI") & " sec.)"
    tmrITI.Interval = datTrialInfo.Recordset("ITI") * 1000
    flagAbort = False
    tmrITI.Enabled = True
End Sub



Private Sub tmrApparatusDelay_Timer()
    Dim Drop As Integer
    
    tmrApparatusDelay.Enabled = False
    tmrApparatusDelay.Interval = 0
    
 If lblApparatus.Caption = "Juice" Then
 
 With DriverLINXSR1
    .Req_subsystem = DL_DO
    .Req_mode = DL_POLLED
    .Req_op = DL_START
    .Evt_Str_type = DL_NULLEVENT
    .Evt_Stp_type = DL_NULLEVENT
    .Evt_Tim_type = DL_NULLEVENT
    .Sel_chan_format = DL_tNATIVE
    .Sel_chan_N = 1
    .Sel_chan_start = 1   'Channel 0 is configured for Output already
    .Sel_chan_startGainCode = 0
    .Sel_buf_samples = 1
    .Sel_buf_N = 0
    .Refresh
    .Res_Sta_ioValue = APPARATUS_JUICE
    .Refresh
    Label1.Caption = .Message
End With

tmrMetraByte.Interval = datSessionInfo.Recordset("JuiceTime")
        tmrMetraByte.Enabled = True
        lblNumJuice.Caption = lblNumJuice.Caption + 1
    
     
              
'        agOutp METRABYTE_PORT_B, APPARATUS_OFF
  '      tmrMetraByte.Interval = datSessionInfo.Recordset("DefaultIncorrectLTO") * 1000 'datSessionInfo.Recordset("JuiceTime")
    '    tmrMetraByte.Interval = 0
    '    tmrMetraByte.Enabled = True
    '    lblNumJuice.Caption = lblNumJuice.Caption + 1
       
    
    ElseIf lblApparatus.Caption = "Pellets" Then
       
   ' If lblApparatus.Caption = "Pellets" Then
       'multiple hopper
        'If datTrialInfo.Recordset("TrialType") = "Left" Then Drop = Int(Rnd * 3)
        'If datTrialInfo.Recordset("TrialType") = "Right" Then Drop = Int(Rnd * 3) + 3
       ' Select Case Drop
       '     Case 0
       '         agOutp METRABYTE_PORT_B, APPARATUS_HOPPER1
       '     Case 1
       '         agOutp METRABYTE_PORT_B, APPARATUS_HOPPER2
       '     Case 2
       '         agOutp METRABYTE_PORT_B, APPARATUS_HOPPER3
       '     Case 3
       '         agOutp METRABYTE_PORT_B, APPARATUS_HOPPER4
       '     Case 4
       '         agOutp METRABYTE_PORT_B, APPARATUS_HOPPER5
       '     Case 5
       '         agOutp METRABYTE_PORT_B, APPARATUS_HOPPER6
       ' End Select
        'agOutp METRABYTE_CONTROL, APPARATUS_SET_B_OUTPUT
        'agOutp METRABYTE_PORT_B, APPARATUS_HOPPER + APPARATUS_HOPPER_LIGHT
        
With DriverLINXSR1
    .Req_subsystem = DL_DO
    .Req_mode = DL_POLLED
    .Req_op = DL_START
    .Evt_Str_type = DL_NULLEVENT
    .Evt_Stp_type = DL_NULLEVENT
    .Evt_Tim_type = DL_NULLEVENT
    .Sel_chan_format = DL_tNATIVE
    .Sel_chan_N = 1
    .Sel_chan_start = 1   'Channel 0 is configured for Output already
    .Sel_chan_startGainCode = 0
    .Sel_buf_samples = 1
    .Sel_buf_N = 0
    .Refresh
    .Res_Sta_ioValue = APPARATUS_PELLETS
    .Refresh
    Label1.Caption = .Message
End With
     
        
        
        
        tmrMetraByte.Interval = datSessionInfo.Recordset("PelletTime")
        tmrMetraByte.Enabled = True
        lblNumPellets.Caption = lblNumPellets.Caption + 1
    
    End If
    
End Sub

Private Sub tmrClock_Timer()

    lblClock.Caption = Format(Now, "hh:mm:ss")
    
End Sub

Private Sub tmrCorrectItemDuration_Timer()
    tmrCorrectItemDuration.Enabled = False
    tmrCorrectItemDuration.Interval = 0
    subClearDisplay
    lblViewTime.Caption = ""
End Sub

Private Sub tmrHouseLight_Timer()

    tmrHouseLight.Enabled = False
    tmrHouseLight.Interval = 0
    
    'agOutp METRABYTE_PORT_B, APPARATUS_OFF
    DriverLINXSR1.Res_Sta_ioValue = APPARATUS_OFF
    DriverLINXSR1.Refresh

End Sub

Private Sub tmrHouseLight1_Timer()

    tmrHouseLight1.Enabled = False
    tmrHouseLight1.Interval = 0
    
    
    
    
    showReadySignal
    
    
    'subInitialItems
    
    listenToFormPecks = True
End Sub

Private Sub tmrIncorrectLTO_Timer()

    tmrIncorrectLTO.Enabled = False
    tmrIncorrectLTO.Interval = 0
    
    lblNumJuice.Caption = lblNumJuice.Caption + 1
'    agOutp METRABYTE_CONTROL, APPARATUS_SET_B_OUTPUT
'    agOutp METRABYTE_PORT_B, APPARATUS_HOUSE_LIGHT

With DriverLINXSR1
    .Req_subsystem = DL_DO
    .Req_mode = DL_POLLED
    .Req_op = DL_START
    .Evt_Str_type = DL_NULLEVENT
    .Evt_Stp_type = DL_NULLEVENT
    .Evt_Tim_type = DL_NULLEVENT
    .Sel_chan_format = DL_tNATIVE
    .Sel_chan_N = 1
    .Sel_chan_start = 1   'Channel 0 is configured for Output already
    .Sel_chan_startGainCode = 0
    .Sel_buf_samples = 1
    .Sel_buf_N = 0
    .Refresh
    .Res_Sta_ioValue = APPARATUS_HOUSE_LIGHT
    .Refresh
    Label1.Caption = .Message
End With

    tmrHouseLight.Interval = (datTrialInfo.Recordset("ITI") * 1000) - 300
    tmrHouseLight.Enabled = True
    
    lblStateMessage.Caption = "ITI (" & datTrialInfo.Recordset("ITI") & " sec.)"
    tmrITI.Interval = datTrialInfo.Recordset("ITI") * 1000
    flagAbort = False
    tmrITI.Enabled = True

    
End Sub


Private Sub tmrInit_Timer()
    
  'mt  subInitTouchscreenRun tmrInit, sInitState, tmrReadTouch

End Sub

Private Sub tmrISI_Timer()
    
    tmrISI.Enabled = False
    tmrISI.Interval = 0
Stop
    subDisplayStimuli
    
End Sub

Private Sub tmrITI_Timer()
    listenToFormPecks = True

    tmrITI.Enabled = False
    tmrITI.Interval = 0
    
    'agOutp METRABYTE_PORT_B, APPARATUS_OFF
    DriverLINXSR1.Res_Sta_ioValue = APPARATUS_OFF
    DriverLINXSR1.Refresh
    
    
    'agOutp METRABYTE_CONTROL, APPARATUS_SET_B_OUTPUT 'initialize
    'initialize PCI driver
    With DriverLINXSR1
        .Req_device = 0
        .Req_mode = DL_OTHER
        .Req_op = DL_INITIALIZE
        .Req_DLL_name = "KPCIPIO"
        .Refresh
        Label1.Caption = .Message
    End With
    
    ' code below will Configure Channel 1 for Output
    With DriverLINXSR1
        .Req_subsystem = DL_DO  ' the subsystem to assign the channel to
        .Req_mode = DL_OTHER
        .Req_op = DL_CONFIGURE  ' it is a configuration type operation
        .Evt_Tim_type = DL_DIOSETUP
        .Evt_Tim_dioChannel = 1   ' configure channel 1 or Port b as Output
        .Evt_Tim_dioMode = DL_DIO_BASIC
        .Evt_Str_type = DL_NULLEVENT
        .Evt_Stp_type = DL_NULLEVENT
        .Sel_chan_N = 0
        .Refresh
        lblStateMessage.Caption = .Message
    End With
    If datSessionInfo.Recordset("IncorrectCorrectionsEnabled") And TrialStatus = False Then
        datStimuliInfo.Recordset.MoveFirst
    '    If flagAbort = False Then
        'subInitialItems
        showReadySignal
    Else
        subTrialResults
        datTrialInfo.Recordset.MoveNext
        If datTrialInfo.Recordset.EOF Then
            subStopSession
        Else
            datStimuliInfo.RecordSource = "SELECT * FROM StimuliInfo WHERE TrialNum = " & datTrialInfo.Recordset("TrialNum") & " ORDER BY Position"
            datStimuliInfo.Refresh
            datProbeInfo.RecordSource = "SELECT * FROM ProbeInfo WHERE TrialNum = " & datTrialInfo.Recordset("TrialNum")
            datProbeInfo.Refresh
            If flagAbort = False Then showReadySignal 'subInitialItems
        End If 'EOF
    End If 'timeout
End Sub

Private Sub tmrMetraByte_Timer()
    
    tmrMetraByte.Enabled = False
    tmrMetraByte.Interval = 0
    
    'agOutp METRABYTE_PORT_B, APPARATUS_OFF
    DriverLINXSR1.Res_Sta_ioValue = APPARATUS_OFF
    DriverLINXSR1.Refresh
    
   ' agOutp METRABYTE_PORT_B, APPARATUS_HOUSE_LIGHT
   With DriverLINXSR1
    .Req_subsystem = DL_DO
    .Req_mode = DL_POLLED
    .Req_op = DL_START
    .Evt_Str_type = DL_NULLEVENT
    .Evt_Stp_type = DL_NULLEVENT
    .Evt_Tim_type = DL_NULLEVENT
    .Sel_chan_format = DL_tNATIVE
    .Sel_chan_N = 1
    .Sel_chan_start = 1   'Channel 0 is configured for Output already
    .Sel_chan_startGainCode = 0
    .Sel_buf_samples = 1
    .Sel_buf_N = 0
    .Refresh
    .Res_Sta_ioValue = APPARATUS_HOUSE_LIGHT
    .Refresh
    Label1.Caption = .Message
End With

    tmrHouseLight.Interval = (datTrialInfo.Recordset("ITI") * 1000) - 300
    tmrHouseLight.Enabled = True

    lblStateMessage.Caption = "ITI (" & datTrialInfo.Recordset("ITI") & " sec.)"
    tmrITI.Interval = datTrialInfo.Recordset("ITI") * 1000
    flagAbort = False
    tmrITI.Enabled = True
    

End Sub


Private Sub tmrProbeAbort_Timer()
    tmrProbeAbort.Enabled = False
    tmrProbeAbort.Interval = 0

    lblStateMessage.Caption = "ITI (" & datTrialInfo.Recordset("ITI") & " sec.)"
    tmrITI.Interval = datTrialInfo.Recordset("ITI") * 1000
    flagAbort = False
    
    'agOutp METRABYTE_PORT_B, APPARATUS_HOUSE_LIGHT
With DriverLINXSR1
    .Req_subsystem = DL_DO
    .Req_mode = DL_POLLED
    .Req_op = DL_START
    .Evt_Str_type = DL_NULLEVENT
    .Evt_Stp_type = DL_NULLEVENT
    .Evt_Tim_type = DL_NULLEVENT
    .Sel_chan_format = DL_tNATIVE
    .Sel_chan_N = 1
    .Sel_chan_start = 1   'Channel 0 is configured for Output already
    .Sel_chan_startGainCode = 0
    .Sel_buf_samples = 1
    .Sel_buf_N = 0
    .Refresh
    .Res_Sta_ioValue = APPARATUS_HOUSE_LIGHT
    .Refresh
    Label1.Caption = .Message
End With
    
    
    tmrITI.Enabled = True

End Sub

Private Sub tmrProbeDelay_Timer()

    tmrProbeDelay.Enabled = False
    tmrProbeDelay.Interval = 0
    
    subSecondItems
    'subDisplayProbe
    
End Sub

Private Sub tmrProbeView_Timer()
    
    lResponseTime = lResponseTime + tmrProbeView.Interval
    lblStateMessage.Caption = lResponseTime & " of " & gDisplayTime
    If lResponseTime >= datSessionInfo.Recordset("ProbeTime") And _
    frmDisplay.leadProbe.Visible = True And datSessionInfo.Recordset("ProbeAbortEnabled") Then sResponse = "Abort"
    
    If datSessionInfo.Recordset("OperationalMode") = "Normal" Then
        subCheckResponse
    End If
        
End Sub

Private Sub tmrReadTouch_Timer()
    
    Dim iX As Integer
    Dim iY As Integer
    Dim adjust As Integer
    Dim adjustd As Integer
    
   ' gsubGetXY iX, iY
    
    adjust = 60
    adjustd = 20
    If lblSubject.Caption = "Aruba" Or lblSubject.Caption = "Gracie" Or lblSubject.Caption = "P2" Or lblSubject.Caption = "Slim" Then
        adjust = 0
        adjustd = 0
    End If
    
    If sImageDisplayed = "stim" Then
    
        If iX >= (STIMULI_LEFT - adjust) And iX <= (STIMULI_LEFT + STIMULI_WIDTH) And iY >= (STIMULI_TOP - 200) And iY <= (STIMULI_TOP + STIMULI_HEIGHT) Then
            sResponse = "Stim"
        Else
            sResponse = "None"
        End If
        
        
    ElseIf sImageDisplayed = "probe" Then
    
        If iX >= (PROBE_LEFT - adjust) And iX <= (PROBE_LEFT + PROBE_WIDTH) And iY >= PROBE_TOP And iY <= (PROBE_TOP + PROBE_HEIGHT) Then
            sResponse = "Same"
            
        ElseIf iX >= (DIFF_LEFT - adjustd) And iX <= ((DIFF_LEFT - adjustd) + DIFF_WIDTH) And iY >= DIFF_TOP And iY <= (DIFF_TOP + DIFF_HEIGHT) Then
            sResponse = "Diff"
    
        Else
            sResponse = "None"
        
        End If
    
    Else
        sResponse = "None"

    End If
    
End Sub

Private Sub tmrSampleDuration_Timer()

    tmrSampleDuration.Enabled = False
    tmrSampleDuration.Interval = 0
    leadStimuli.Visible = False
    frmDisplay.leadStimuli.Visible = False
    
End Sub

Private Sub tmrStimuliView_Timer()
    
    Dim lISITime As Long
    
    tmrStimuliView.Enabled = False
    tmrStimuliView.Interval = 0
 
    If datSessionInfo.Recordset("Phase") = 2 Or datSessionInfo.Recordset("Phase") = 3 Then
        leadStimuli.Visible = False
        frmDisplay.leadStimuli.Visible = False
    End If
    frmDisplay.leadProbe.Visible = False
   ' frmDisplay.leadProbe2.Visible = False
    frmDisplay.lblDiff.Visible = False
    'subClearDisplay
     
     
'If iStimTouches = 0 And datSessionInfo.Recordset("TouchAbortEnabled") And datStimuliInfo.Recordset("Position") = 1 Then    'only aborts on first list item
'    subResetProbe
'    subAbortLTO
'Else
    lISITime = datStimuliInfo.Recordset("ISI")
    datStimuliInfo.Recordset.MoveNext
    
    If datStimuliInfo.Recordset.EOF Then
        tmrISI.Enabled = False
        tmrISI.Interval = 0
        
        lblStateMessage.Caption = "Probe Delay (" & datProbeInfo.Recordset("Delay") & " sec.)"
        If datProbeInfo.Recordset("Delay") = 0 Then
            subDisplayProbe
        Else
            tmrProbeDelay.Interval = (datProbeInfo.Recordset("Delay") * 1000)
            tmrProbeDelay.Enabled = True
        End If
    Else
        lblStateMessage.Caption = "ISI (" & lISITime & " sec.)"
        tmrISI.Interval = lISITime * 1000
        tmrISI.Enabled = True

    End If
'End If 'abort
End Sub



Public Sub subRecordItem()

Select Case datStimuliInfo.Recordset("Position")
    Case 1
        Item1 = datStimuliInfo.Recordset("Filename")
        Item2 = "blank"
        Item3 = "blank"
        Item4 = "blank"
        Item5 = "blank"
        Item6 = "blank"
        Item7 = "blank"
        Item8 = "blank"
        Item9 = "blank"
        Item10 = "blank"
    Case 2
        Item2 = datStimuliInfo.Recordset("Filename")
    Case 3
        Item3 = datStimuliInfo.Recordset("Filename")
    Case 4
        Item4 = datStimuliInfo.Recordset("Filename")
    Case 5
        Item5 = datStimuliInfo.Recordset("Filename")
    Case 6
        Item6 = datStimuliInfo.Recordset("Filename")
    Case 7
        Item7 = datStimuliInfo.Recordset("Filename")
    Case 8
        Item8 = datStimuliInfo.Recordset("Filename")
    Case 9
        Item9 = datStimuliInfo.Recordset("Filename")
    Case 10
        Item10 = datStimuliInfo.Recordset("Filename")
End Select
End Sub

Public Sub subRecordTouch()

    gPecks = gPecks + 1
    datStimResponses.Recordset.AddNew
    datStimResponses.Recordset("TrialNum") = lblTrialNum.Caption
    datStimResponses.Recordset("CorrectionTrial") = lblTotalIncorrectCorrections.Caption
    datStimResponses.Recordset("ListPosition") = RespLoc
    datStimResponses.Recordset("ResponseTime") = lResponseTime
    datStimResponses.Recordset.Update

End Sub
