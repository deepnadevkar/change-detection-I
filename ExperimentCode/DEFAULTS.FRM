VERSION 5.00
Object = "{02B97C13-ED1E-11CD-A08B-00AA00575482}#1.0#0"; "mhiinp32.ocx"
Object = "{B08ECCE3-EE86-11CE-9577-0020AF039CA3}#1.0#0"; "mhwave32.ocx"
Object = "{F9043C88-F6F2-101A-A3C9-08002B2F49FB}#1.2#0"; "COMDLG32.OCX"
Begin VB.Form frmDefaults 
   BackColor       =   &H80000000&
   BorderStyle     =   1  'Fixed Single
   Caption         =   "Change Detection Defaults"
   ClientHeight    =   8010
   ClientLeft      =   225
   ClientTop       =   510
   ClientWidth     =   10050
   ControlBox      =   0   'False
   LinkTopic       =   "Form1"
   MaxButton       =   0   'False
   MinButton       =   0   'False
   PaletteMode     =   1  'UseZOrder
   ScaleHeight     =   8010
   ScaleWidth      =   10050
   Begin MhinintLibCtl.MhIntInput IntProbeAbortLTO 
      DataField       =   "ProbeAbortLTO"
      DataSource      =   "datDefaults"
      Height          =   285
      Left            =   1440
      TabIndex        =   104
      Top             =   3600
      Width           =   465
      _Version        =   65536
      _ExtentX        =   820
      _ExtentY        =   503
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput IntProbeTime 
      DataField       =   "ProbeTime"
      DataSource      =   "datDefaults"
      Height          =   285
      Left            =   840
      TabIndex        =   99
      Top             =   3600
      Width           =   465
      _Version        =   65536
      _ExtentX        =   820
      _ExtentY        =   503
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput VT 
      DataField       =   "pd8"
      DataSource      =   "datView"
      Height          =   315
      Index           =   10
      Left            =   4680
      TabIndex        =   48
      Top             =   3120
      Width           =   555
      _Version        =   65536
      _ExtentX        =   979
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   16777215
      BevelSize       =   2
      BorderStyle     =   0
      FillColor       =   16777215
      ShadowColor     =   -2147483630
      Min             =   -2147483647
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput VT 
      DataField       =   "pd4"
      DataSource      =   "datView"
      Height          =   315
      Index           =   9
      Left            =   4680
      TabIndex        =   44
      Top             =   2760
      Width           =   555
      _Version        =   65536
      _ExtentX        =   979
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   16777215
      BevelSize       =   2
      BorderStyle     =   0
      FillColor       =   16777215
      ShadowColor     =   -2147483630
      Min             =   -2147483647
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput VT 
      DataField       =   "pd7"
      DataSource      =   "datView"
      Height          =   315
      Index           =   8
      Left            =   4080
      TabIndex        =   47
      Top             =   3120
      Width           =   555
      _Version        =   65536
      _ExtentX        =   979
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   16777215
      BevelSize       =   2
      BorderStyle     =   0
      FillColor       =   16777215
      ShadowColor     =   -2147483630
      Min             =   -2147483647
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput VT 
      DataField       =   "pd3"
      DataSource      =   "datView"
      Height          =   315
      Index           =   7
      Left            =   4080
      TabIndex        =   43
      Top             =   2760
      Width           =   555
      _Version        =   65536
      _ExtentX        =   979
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   16777215
      BevelSize       =   2
      BorderStyle     =   0
      FillColor       =   16777215
      ShadowColor     =   -2147483630
      Min             =   -2147483647
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput VT 
      DataField       =   "ds7"
      DataSource      =   "datView"
      Height          =   285
      Index           =   6
      Left            =   4320
      TabIndex        =   55
      Top             =   5400
      Width           =   615
      _Version        =   65536
      _ExtentX        =   1085
      _ExtentY        =   503
      _StockProps     =   77
      BackColor       =   16777215
      BevelSize       =   2
      BorderStyle     =   0
      FillColor       =   16777215
      Min             =   -2147483647
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput VT 
      DataField       =   "ds6"
      DataSource      =   "datView"
      Height          =   285
      Index           =   5
      Left            =   4320
      TabIndex        =   54
      Top             =   5040
      Width           =   615
      _Version        =   65536
      _ExtentX        =   1085
      _ExtentY        =   503
      _StockProps     =   77
      BackColor       =   16777215
      BevelSize       =   2
      BorderStyle     =   0
      FillColor       =   16777215
      Min             =   -2147483647
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput VT 
      DataField       =   "ds5"
      DataSource      =   "datView"
      Height          =   285
      Index           =   4
      Left            =   3240
      TabIndex        =   53
      Top             =   6000
      Width           =   615
      _Version        =   65536
      _ExtentX        =   1085
      _ExtentY        =   503
      _StockProps     =   77
      BackColor       =   16777215
      BevelSize       =   2
      BorderStyle     =   0
      FillColor       =   16777215
      Min             =   -2147483647
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput VT 
      DataField       =   "ds4"
      DataSource      =   "datView"
      Height          =   285
      Index           =   3
      Left            =   3240
      TabIndex        =   52
      Top             =   5640
      Width           =   615
      _Version        =   65536
      _ExtentX        =   1085
      _ExtentY        =   503
      _StockProps     =   77
      BackColor       =   16777215
      BevelSize       =   2
      BorderStyle     =   0
      FillColor       =   16777215
      Min             =   -2147483647
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput VT 
      DataField       =   "ds3"
      DataSource      =   "datView"
      Height          =   285
      Index           =   2
      Left            =   3240
      TabIndex        =   51
      Top             =   5340
      Width           =   615
      _Version        =   65536
      _ExtentX        =   1085
      _ExtentY        =   503
      _StockProps     =   77
      BackColor       =   16777215
      BevelSize       =   2
      BorderStyle     =   0
      FillColor       =   16777215
      Min             =   -2147483647
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput VT 
      DataField       =   "ds2"
      DataSource      =   "datView"
      Height          =   285
      Index           =   1
      Left            =   3240
      TabIndex        =   50
      Top             =   5040
      Width           =   615
      _Version        =   65536
      _ExtentX        =   1085
      _ExtentY        =   503
      _StockProps     =   77
      BackColor       =   16777215
      BevelSize       =   2
      BorderStyle     =   0
      FillColor       =   16777215
      Min             =   -2147483647
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intStartTrial 
      DataField       =   "StartTrial"
      DataSource      =   "datDefaults"
      Height          =   315
      Left            =   1620
      TabIndex        =   6
      Top             =   1920
      Width           =   1065
      _Version        =   65536
      _ExtentX        =   1879
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Text            =   "1"
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intProbeTouches 
      DataField       =   "ProbeTouches"
      DataSource      =   "datDefaults"
      Height          =   315
      Index           =   1
      Left            =   6360
      TabIndex        =   80
      Top             =   7650
      Visible         =   0   'False
      Width           =   705
      _Version        =   65536
      _ExtentX        =   1244
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Text            =   "1"
      Max             =   20
      Min             =   1
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intStimTouches 
      DataField       =   "StimTouches"
      DataSource      =   "datDefaults"
      Height          =   315
      Index           =   1
      Left            =   1620
      TabIndex        =   11
      Top             =   4320
      Width           =   705
      _Version        =   65536
      _ExtentX        =   1244
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Text            =   "1"
      Max             =   20
      Min             =   1
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhwaveLib.MhWave waveTest 
      Height          =   420
      Left            =   7950
      TabIndex        =   77
      Top             =   360
      Visible         =   0   'False
      Width           =   420
      _Version        =   65536
      _ExtentX        =   741
      _ExtentY        =   741
      _StockProps     =   0
      FilenameRead    =   "DEFAULTS.frx":0000
   End
   Begin MhinintLibCtl.MhIntInput intJuiceTime 
      DataField       =   "JuiceTime"
      DataSource      =   "datDefaults"
      Height          =   315
      Index           =   0
      Left            =   1620
      TabIndex        =   70
      Top             =   6390
      Visible         =   0   'False
      Width           =   705
      _Version        =   65536
      _ExtentX        =   1244
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
   Begin MhinintLibCtl.MhIntInput intPelletTime 
      DataField       =   "PelletTime"
      DataSource      =   "datDefaults"
      Height          =   315
      Index           =   0
      Left            =   1620
      TabIndex        =   13
      Top             =   5280
      Width           =   705
      _Version        =   65536
      _ExtentX        =   1244
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Text            =   "2000"
      Max             =   30000
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intDiffProb 
      DataField       =   "pd6"
      DataSource      =   "datView"
      Height          =   315
      Left            =   3480
      TabIndex        =   46
      Top             =   3120
      Width           =   555
      _Version        =   65536
      _ExtentX        =   979
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Max             =   100
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intSameProb 
      DataField       =   "pd5"
      DataSource      =   "datView"
      Height          =   315
      Left            =   2880
      TabIndex        =   45
      Top             =   3120
      Width           =   555
      _Version        =   65536
      _ExtentX        =   979
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Max             =   100
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intReinfProb 
      DataField       =   "ReinfProbability"
      DataSource      =   "datDefaults"
      Height          =   315
      Left            =   8040
      TabIndex        =   42
      Top             =   2520
      Width           =   555
      _Version        =   65536
      _ExtentX        =   979
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Max             =   100
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intPelletProb 
      DataField       =   "PelletProbability"
      DataSource      =   "datDefaults"
      Height          =   315
      Left            =   4680
      TabIndex        =   41
      Top             =   2160
      Width           =   555
      _Version        =   65536
      _ExtentX        =   979
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Text            =   "100"
      Max             =   100
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intJuiceProb 
      DataField       =   "JuiceProbability"
      DataSource      =   "datDefaults"
      Height          =   315
      Left            =   1650
      TabIndex        =   39
      Top             =   6780
      Visible         =   0   'False
      Width           =   555
      _Version        =   65536
      _ExtentX        =   979
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Max             =   100
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intNumTrials 
      DataField       =   "NumOfTrials"
      DataSource      =   "datDefaults"
      Height          =   315
      Left            =   2340
      TabIndex        =   5
      Top             =   1470
      Width           =   345
      _Version        =   65536
      _ExtentX        =   609
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Text            =   "96"
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intListLength 
      DataField       =   "DefaultMaxListLength"
      DataSource      =   "datDefaults"
      Height          =   315
      Left            =   9240
      TabIndex        =   36
      Top             =   6960
      Visible         =   0   'False
      Width           =   705
      _Version        =   65536
      _ExtentX        =   1244
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Text            =   "1"
      Min             =   1
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intViewTime 
      DataField       =   "DefaultViewTime"
      DataSource      =   "datDefaults"
      Height          =   315
      Left            =   6330
      TabIndex        =   34
      Top             =   7320
      Visible         =   0   'False
      Width           =   705
      _Version        =   65536
      _ExtentX        =   1244
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intITI 
      DataField       =   "DefaultITI"
      DataSource      =   "datDefaults"
      Height          =   285
      Left            =   720
      TabIndex        =   9
      Top             =   3960
      Visible         =   0   'False
      Width           =   465
      _Version        =   65536
      _ExtentX        =   820
      _ExtentY        =   503
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intISI 
      DataField       =   "DefaultISI"
      DataSource      =   "datDefaults"
      Height          =   315
      Left            =   1560
      TabIndex        =   31
      Top             =   5760
      Width           =   705
      _Version        =   65536
      _ExtentX        =   1244
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intProbeDelay 
      DataField       =   "DefaultProbeDelay"
      DataSource      =   "datDefaults"
      Height          =   315
      Left            =   1620
      TabIndex        =   10
      Top             =   3900
      Visible         =   0   'False
      Width           =   705
      _Version        =   65536
      _ExtentX        =   1244
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Text            =   "50"
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intAbortLTO 
      DataField       =   "DefaultAbortLTO"
      DataSource      =   "datDefaults"
      Height          =   285
      Left            =   1500
      TabIndex        =   17
      Top             =   7650
      Visible         =   0   'False
      Width           =   705
      _Version        =   65536
      _ExtentX        =   1244
      _ExtentY        =   503
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intIncorrectLTO 
      DataField       =   "DefaultIncorrectLTO"
      DataSource      =   "datDefaults"
      Height          =   315
      Left            =   1620
      TabIndex        =   8
      Top             =   2730
      Width           =   705
      _Version        =   65536
      _ExtentX        =   1244
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput VT 
      DataField       =   "pd2"
      DataSource      =   "datView"
      Height          =   315
      Index           =   0
      Left            =   3480
      TabIndex        =   113
      Top             =   2760
      Width           =   555
      _Version        =   65536
      _ExtentX        =   979
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   16777215
      BevelSize       =   2
      BorderStyle     =   0
      FillColor       =   16777215
      ShadowColor     =   -2147483630
      Min             =   -2147483647
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput VT 
      DataField       =   "pd1"
      DataSource      =   "datView"
      Height          =   315
      Index           =   11
      Left            =   2880
      TabIndex        =   114
      Top             =   2760
      Width           =   555
      _Version        =   65536
      _ExtentX        =   979
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   16777215
      BevelSize       =   2
      BorderStyle     =   0
      FillColor       =   16777215
      ShadowColor     =   -2147483630
      Min             =   -2147483647
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput VT 
      DataField       =   "ds8"
      DataSource      =   "datView"
      Height          =   285
      Index           =   12
      Left            =   4320
      TabIndex        =   116
      Top             =   5760
      Width           =   615
      _Version        =   65536
      _ExtentX        =   1085
      _ExtentY        =   503
      _StockProps     =   77
      BackColor       =   16777215
      BevelSize       =   2
      BorderStyle     =   0
      FillColor       =   16777215
      Min             =   -2147483647
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput VT 
      DataField       =   "ds9"
      DataSource      =   "datView"
      Height          =   285
      Index           =   13
      Left            =   4320
      TabIndex        =   118
      Top             =   6060
      Width           =   615
      _Version        =   65536
      _ExtentX        =   1085
      _ExtentY        =   503
      _StockProps     =   77
      BackColor       =   16777215
      BevelSize       =   2
      BorderStyle     =   0
      FillColor       =   16777215
      Min             =   -2147483647
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput VT 
      DataField       =   "ds10"
      DataSource      =   "datView"
      Height          =   285
      Index           =   14
      Left            =   5640
      TabIndex        =   120
      Top             =   5040
      Width           =   615
      _Version        =   65536
      _ExtentX        =   1085
      _ExtentY        =   503
      _StockProps     =   77
      BackColor       =   16777215
      BevelSize       =   2
      BorderStyle     =   0
      FillColor       =   16777215
      Min             =   -2147483647
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput VT 
      DataField       =   "ds11"
      DataSource      =   "datView"
      Height          =   285
      Index           =   15
      Left            =   5640
      TabIndex        =   122
      Top             =   5400
      Width           =   615
      _Version        =   65536
      _ExtentX        =   1085
      _ExtentY        =   503
      _StockProps     =   77
      BackColor       =   16777215
      BevelSize       =   2
      BorderStyle     =   0
      FillColor       =   16777215
      Min             =   -2147483647
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput VT 
      DataField       =   "ds12"
      DataSource      =   "datView"
      Height          =   285
      Index           =   16
      Left            =   5640
      TabIndex        =   124
      Top             =   5760
      Width           =   615
      _Version        =   65536
      _ExtentX        =   1085
      _ExtentY        =   503
      _StockProps     =   77
      BackColor       =   16777215
      BevelSize       =   2
      BorderStyle     =   0
      FillColor       =   16777215
      Min             =   -2147483647
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin VB.ComboBox Combo3 
      DataField       =   "SesType"
      DataSource      =   "datDefaults"
      Height          =   315
      ItemData        =   "DEFAULTS.frx":B820
      Left            =   1560
      List            =   "DEFAULTS.frx":B82D
      TabIndex        =   110
      Top             =   6120
      Visible         =   0   'False
      Width           =   1065
   End
   Begin VB.CheckBox chkTrialUnique 
      Alignment       =   1  'Right Justify
      Caption         =   "Check for Blocked DS (Ascending seq)"
      DataField       =   "TrialUnique"
      DataSource      =   "datDefaults"
      Height          =   255
      Left            =   2880
      TabIndex        =   107
      Top             =   6360
      Width           =   3135
   End
   Begin VB.CheckBox chkVariableProbeDelay 
      Caption         =   "Check for Blocked Intervals (uncheck for Random) Req. Rebuild"
      DataField       =   "VariableProbeDelay"
      DataSource      =   "datDefaults"
      Height          =   495
      Left            =   2880
      TabIndex        =   49
      Top             =   3480
      Width           =   3135
   End
   Begin VB.CheckBox chkProbeAbort 
      DataField       =   "ProbeAbortEnabled"
      DataSource      =   "datDefaults"
      Height          =   255
      Left            =   4530
      TabIndex        =   101
      Top             =   7140
      Visible         =   0   'False
      Width           =   285
   End
   Begin VB.CheckBox chkNTAbort 
      DataField       =   "TouchAbortEnabled"
      DataSource      =   "datDefaults"
      Height          =   225
      Left            =   2010
      TabIndex        =   96
      Top             =   7410
      Visible         =   0   'False
      Width           =   285
   End
   Begin VB.Data datView 
      Caption         =   "View"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   6630
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   6270
      Visible         =   0   'False
      Width           =   2295
   End
   Begin VB.Data datDefaults 
      Caption         =   "Defaults"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   6600
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   5880
      Visible         =   0   'False
      Width           =   2355
   End
   Begin VB.CheckBox chkNoTouch 
      DataField       =   "NoTouchProbes"
      DataSource      =   "datDefaults"
      Height          =   255
      Left            =   2130
      TabIndex        =   88
      Top             =   4620
      Visible         =   0   'False
      Width           =   255
   End
   Begin VB.ComboBox Combo2 
      DataField       =   "SetSize"
      DataSource      =   "datDefaults"
      Height          =   315
      ItemData        =   "DEFAULTS.frx":B845
      Left            =   1620
      List            =   "DEFAULTS.frx":B84F
      TabIndex        =   7
      Top             =   2340
      Width           =   1065
   End
   Begin VB.ComboBox Combo1 
      DataField       =   "Subject"
      DataSource      =   "datDefaults"
      Height          =   315
      ItemData        =   "DEFAULTS.frx":B859
      Left            =   1110
      List            =   "DEFAULTS.frx":B869
      TabIndex        =   2
      Top             =   120
      Width           =   1575
   End
   Begin VB.TextBox txtSessionFile 
      Height          =   285
      Left            =   1110
      TabIndex        =   4
      Top             =   1050
      Width           =   1575
   End
   Begin VB.OptionButton optPhase 
      Caption         =   "Phase 3"
      Height          =   315
      Index           =   2
      Left            =   6810
      TabIndex        =   83
      Top             =   900
      Visible         =   0   'False
      Width           =   945
   End
   Begin VB.OptionButton optPhase 
      Caption         =   "Phase 2"
      Height          =   315
      Index           =   1
      Left            =   5580
      TabIndex        =   82
      Top             =   900
      Visible         =   0   'False
      Width           =   945
   End
   Begin VB.OptionButton optPhase 
      Caption         =   "Phase 1"
      Height          =   315
      Index           =   0
      Left            =   4350
      TabIndex        =   81
      Top             =   900
      Width           =   945
   End
   Begin VB.CommandButton cmdTestSound 
      Caption         =   "Test"
      Enabled         =   0   'False
      Height          =   315
      Left            =   8430
      TabIndex        =   76
      Top             =   1350
      Width           =   585
   End
   Begin VB.CommandButton cmdSoundFile 
      Caption         =   "Sound File"
      Height          =   315
      Left            =   3000
      TabIndex        =   75
      Top             =   1350
      Width           =   1125
   End
   Begin VB.CommandButton cmdTestMode 
      Caption         =   "Test - 1 Monitor"
      Height          =   315
      Left            =   4350
      TabIndex        =   63
      Top             =   480
      Width           =   1425
   End
   Begin VB.CommandButton cmdNormalMode 
      Caption         =   "Normal - 2 Monitors"
      Height          =   315
      Left            =   5790
      TabIndex        =   62
      Top             =   480
      Width           =   1695
   End
   Begin VB.CheckBox Check1 
      DataField       =   "DefaultVariableListLength"
      DataSource      =   "datDefaults"
      Enabled         =   0   'False
      Height          =   225
      Left            =   8790
      TabIndex        =   60
      Top             =   7650
      Visible         =   0   'False
      Width           =   315
   End
   Begin VB.DirListBox dirImages 
      Height          =   2790
      Left            =   6540
      TabIndex        =   56
      Top             =   2880
      Width           =   2475
   End
   Begin VB.TextBox txtSessionID 
      DataField       =   "SessionID"
      DataSource      =   "datDefaults"
      Height          =   285
      Left            =   1110
      TabIndex        =   3
      Text            =   "1"
      Top             =   600
      Width           =   1575
   End
   Begin VB.CheckBox chkAbortCorrections 
      DataField       =   "AbortCorrectionsEnabled"
      DataSource      =   "datDefaults"
      Height          =   255
      Left            =   9000
      TabIndex        =   14
      Top             =   7290
      Visible         =   0   'False
      Width           =   285
   End
   Begin VB.CheckBox chkIncorrectCorrections 
      DataField       =   "IncorrectCorrectionsEnabled"
      DataSource      =   "datDefaults"
      Height          =   225
      Left            =   2130
      TabIndex        =   12
      Top             =   4890
      Width           =   285
   End
   Begin VB.CommandButton cmdCancel 
      Caption         =   "Cancel"
      Height          =   345
      Left            =   5040
      TabIndex        =   0
      Top             =   6720
      Width           =   1245
   End
   Begin MSComDlg.CommonDialog diaSound 
      Left            =   7980
      Top             =   810
      _ExtentX        =   847
      _ExtentY        =   847
      _Version        =   393216
   End
   Begin VB.CommandButton cmdOK 
      Caption         =   "OK"
      Height          =   345
      Left            =   3120
      TabIndex        =   1
      Top             =   6720
      Width           =   1245
   End
   Begin VB.Label VTL 
      Alignment       =   1  'Right Justify
      Caption         =   "12"
      Height          =   285
      Index           =   10
      Left            =   5160
      TabIndex        =   125
      Top             =   5760
      Width           =   435
   End
   Begin VB.Label VTL 
      Alignment       =   1  'Right Justify
      Caption         =   "11"
      Height          =   285
      Index           =   9
      Left            =   5280
      TabIndex        =   123
      Top             =   5400
      Width           =   315
   End
   Begin VB.Label VTL 
      Alignment       =   1  'Right Justify
      Caption         =   "10"
      Height          =   285
      Index           =   8
      Left            =   5280
      TabIndex        =   121
      Top             =   5040
      Width           =   315
   End
   Begin VB.Label VTL 
      Alignment       =   1  'Right Justify
      Caption         =   "9"
      Height          =   285
      Index           =   7
      Left            =   3960
      TabIndex        =   119
      Top             =   6060
      Width           =   315
   End
   Begin VB.Label VTL 
      Alignment       =   1  'Right Justify
      Caption         =   "8"
      Height          =   285
      Index           =   0
      Left            =   3960
      TabIndex        =   117
      Top             =   5760
      Width           =   315
   End
   Begin VB.Line Line4 
      X1              =   2760
      X2              =   9060
      Y1              =   2160
      Y2              =   2160
   End
   Begin VB.Label Label25 
      Caption         =   "Juice is 100 - Pellet"
      Height          =   375
      Left            =   5280
      TabIndex        =   115
      Top             =   2160
      Width           =   1455
   End
   Begin VB.Label Label24 
      Caption         =   "Pellet Probability (0-100)"
      Height          =   255
      Left            =   2880
      TabIndex        =   112
      Top             =   2220
      Width           =   1815
   End
   Begin VB.Line Line6 
      X1              =   2790
      X2              =   6360
      Y1              =   4800
      Y2              =   4800
   End
   Begin VB.Label Label23 
      Alignment       =   2  'Center
      BackColor       =   &H80000005&
      Caption         =   "All Retention intervals boxes MUST have a nonnegative value. Use duplicate values to change probability of a given interval"
      Height          =   735
      Left            =   2880
      TabIndex        =   111
      Top             =   3960
      Width           =   3135
   End
   Begin VB.Label Label51 
      Caption         =   "Session Type:"
      Height          =   255
      Left            =   120
      TabIndex        =   109
      Top             =   6150
      Visible         =   0   'False
      Width           =   1095
   End
   Begin VB.Line Line7 
      X1              =   -240
      X2              =   9120
      Y1              =   6840
      Y2              =   6840
   End
   Begin VB.Line Line5 
      X1              =   9120
      X2              =   9120
      Y1              =   0
      Y2              =   6720
   End
   Begin VB.Label Label50 
      Caption         =   "Trial Unique:"
      Height          =   225
      Left            =   6960
      TabIndex        =   108
      Top             =   7050
      Visible         =   0   'False
      Width           =   1605
   End
   Begin VB.Label Label49 
      Caption         =   "Variable Probe Delays:"
      Height          =   225
      Left            =   3720
      TabIndex        =   106
      Top             =   6840
      Visible         =   0   'False
      Width           =   1605
   End
   Begin VB.Label Label48 
      Alignment       =   1  'Right Justify
      Caption         =   "ProbeAbortLTO:"
      Height          =   255
      Left            =   2520
      TabIndex        =   105
      Top             =   7740
      Visible         =   0   'False
      Width           =   1155
   End
   Begin VB.Label Label47 
      Caption         =   "msec."
      Height          =   255
      Left            =   4740
      TabIndex        =   103
      Top             =   7710
      Visible         =   0   'False
      Width           =   405
   End
   Begin VB.Label Label46 
      Alignment       =   1  'Right Justify
      Caption         =   "Abort Correction Enabled:"
      Height          =   225
      Left            =   6480
      TabIndex        =   102
      Top             =   7320
      Visible         =   0   'False
      Width           =   2475
   End
   Begin VB.Label Label45 
      Caption         =   "msec."
      Height          =   255
      Left            =   4740
      TabIndex        =   100
      Top             =   7470
      Visible         =   0   'False
      Width           =   405
   End
   Begin VB.Label Label44 
      Alignment       =   1  'Right Justify
      Caption         =   "Probe Time:"
      Height          =   255
      Left            =   2520
      TabIndex        =   98
      Top             =   7440
      Visible         =   0   'False
      Width           =   855
   End
   Begin VB.Label Label43 
      Alignment       =   1  'Right Justify
      Caption         =   "Touch Aborts Enabled:"
      Height          =   225
      Left            =   0
      TabIndex        =   97
      Top             =   7410
      Visible         =   0   'False
      Width           =   1635
   End
   Begin VB.Label VTL 
      Alignment       =   1  'Right Justify
      Caption         =   "7"
      Height          =   285
      Index           =   6
      Left            =   3960
      TabIndex        =   95
      Top             =   5400
      Width           =   315
   End
   Begin VB.Label VTL 
      Alignment       =   1  'Right Justify
      Caption         =   "6"
      Height          =   285
      Index           =   5
      Left            =   3960
      TabIndex        =   94
      Top             =   5040
      Width           =   315
   End
   Begin VB.Label VTL 
      Alignment       =   1  'Right Justify
      Caption         =   "5"
      Height          =   285
      Index           =   4
      Left            =   2880
      TabIndex        =   93
      Top             =   6000
      Width           =   315
   End
   Begin VB.Label VTL 
      Alignment       =   1  'Right Justify
      Caption         =   "4"
      Height          =   285
      Index           =   3
      Left            =   2880
      TabIndex        =   92
      Top             =   5640
      Width           =   315
   End
   Begin VB.Label VTL 
      Alignment       =   1  'Right Justify
      Caption         =   "3"
      Height          =   285
      Index           =   2
      Left            =   2880
      TabIndex        =   91
      Top             =   5340
      Width           =   315
   End
   Begin VB.Label VTL 
      Alignment       =   1  'Right Justify
      Caption         =   "2"
      Height          =   285
      Index           =   1
      Left            =   2880
      TabIndex        =   90
      Top             =   5040
      Width           =   315
   End
   Begin VB.Label Label42 
      Alignment       =   2  'Center
      BackColor       =   &H80000004&
      Caption         =   "Must add to Num Trials"
      Height          =   225
      Left            =   3120
      TabIndex        =   89
      Top             =   4800
      Width           =   2895
   End
   Begin VB.Label Label41 
      Caption         =   "FR-0 Session:"
      Height          =   255
      Left            =   120
      TabIndex        =   87
      Top             =   4680
      Visible         =   0   'False
      Width           =   1455
   End
   Begin VB.Label Label40 
      Caption         =   "Set Size:"
      Height          =   255
      Left            =   120
      TabIndex        =   86
      Top             =   2400
      Width           =   735
   End
   Begin VB.Label Label39 
      Caption         =   "Start at Trial:"
      Height          =   375
      Left            =   120
      TabIndex        =   85
      Top             =   1920
      Width           =   1215
   End
   Begin VB.Label Label38 
      Caption         =   "Session File:"
      Height          =   255
      Left            =   120
      TabIndex        =   84
      Top             =   1080
      Width           =   975
   End
   Begin VB.Label Label36 
      BackColor       =   &H00C0C0C0&
      Caption         =   "Probe Touches:"
      Height          =   255
      Index           =   1
      Left            =   5040
      TabIndex        =   79
      Top             =   7680
      Visible         =   0   'False
      Width           =   1215
   End
   Begin VB.Label Label37 
      BackColor       =   &H80000000&
      Caption         =   "Fixed Ratio:"
      Height          =   255
      Index           =   1
      Left            =   120
      TabIndex        =   78
      Top             =   4320
      Width           =   1095
   End
   Begin VB.Line Line3 
      X1              =   6480
      X2              =   6480
      Y1              =   2520
      Y2              =   6360
   End
   Begin VB.Label lblSound 
      BorderStyle     =   1  'Fixed Single
      DataField       =   "SoundFile"
      DataSource      =   "datDefaults"
      Height          =   315
      Left            =   4350
      TabIndex        =   74
      Top             =   1350
      Width           =   4035
   End
   Begin VB.Label lblPhase 
      Alignment       =   2  'Center
      BorderStyle     =   1  'Fixed Single
      Caption         =   "1"
      DataField       =   "Phase"
      DataSource      =   "datDefaults"
      BeginProperty Font 
         Name            =   "Arial"
         Size            =   12
         Charset         =   0
         Weight          =   700
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Height          =   375
      Left            =   3000
      TabIndex        =   73
      Top             =   840
      Width           =   1215
   End
   Begin VB.Label Label37 
      BackColor       =   &H00C0C0C0&
      Caption         =   "Juice Time:"
      Height          =   255
      Index           =   0
      Left            =   120
      TabIndex        =   72
      Top             =   5760
      Width           =   825
   End
   Begin VB.Label Label36 
      BackColor       =   &H80000000&
      Caption         =   "Grain Access:"
      Height          =   255
      Index           =   0
      Left            =   120
      TabIndex        =   71
      Top             =   5280
      Width           =   1095
   End
   Begin VB.Label Label35 
      Caption         =   "msec."
      Height          =   255
      Left            =   6120
      TabIndex        =   69
      Top             =   7080
      Visible         =   0   'False
      Width           =   405
   End
   Begin VB.Label Label34 
      Caption         =   "msec."
      Height          =   255
      Left            =   2340
      TabIndex        =   68
      Top             =   5310
      Width           =   405
   End
   Begin VB.Label Label33 
      Caption         =   "Results File:"
      Height          =   255
      Left            =   3000
      TabIndex        =   67
      Top             =   1830
      Width           =   915
   End
   Begin VB.Label lblResultsFile 
      BorderStyle     =   1  'Fixed Single
      Height          =   315
      Left            =   4350
      TabIndex        =   66
      Top             =   1770
      Width           =   4665
   End
   Begin VB.Label Label32 
      Caption         =   "Operational Mode:"
      Height          =   375
      Left            =   3030
      TabIndex        =   65
      Top             =   60
      Width           =   1065
   End
   Begin VB.Label lblOpsMode 
      Alignment       =   2  'Center
      BorderStyle     =   1  'Fixed Single
      Caption         =   "Normal"
      DataField       =   "OperationalMode"
      DataSource      =   "datDefaults"
      BeginProperty Font 
         Name            =   "Arial"
         Size            =   12
         Charset         =   0
         Weight          =   700
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Height          =   315
      Left            =   3000
      TabIndex        =   64
      Top             =   480
      Width           =   1215
   End
   Begin VB.Label Label31 
      Caption         =   "Variable ListLength:"
      Height          =   225
      Left            =   7170
      TabIndex        =   61
      Top             =   7680
      Visible         =   0   'False
      Width           =   1515
   End
   Begin VB.Label lblDate 
      BorderStyle     =   1  'Fixed Single
      DataField       =   "Date"
      DataSource      =   "datDefaults"
      Height          =   285
      Left            =   5580
      TabIndex        =   59
      Top             =   60
      Width           =   1575
   End
   Begin VB.Label Label1 
      Caption         =   "Date:"
      Height          =   195
      Left            =   5160
      TabIndex        =   58
      Top             =   90
      Width           =   375
   End
   Begin VB.Label Label14 
      Caption         =   "Image Path:"
      Height          =   195
      Left            =   6450
      TabIndex        =   57
      Top             =   2640
      Width           =   975
   End
   Begin VB.Line Line2 
      X1              =   2880
      X2              =   9150
      Y1              =   6360
      Y2              =   6360
   End
   Begin VB.Line Line1 
      X1              =   2790
      X2              =   2790
      Y1              =   120
      Y2              =   6840
   End
   Begin VB.Label Label26 
      Caption         =   "%"
      Height          =   255
      Left            =   2250
      TabIndex        =   40
      Top             =   6840
      Visible         =   0   'False
      Width           =   345
   End
   Begin VB.Label Label22 
      BackColor       =   &H80000000&
      Caption         =   "Retention Intervals"
      BeginProperty Font 
         Name            =   "MS Sans Serif"
         Size            =   8.25
         Charset         =   0
         Weight          =   700
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Height          =   255
      Left            =   2880
      TabIndex        =   38
      Top             =   2520
      Width           =   1905
   End
   Begin VB.Label Label21 
      BackColor       =   &H00C0C0C0&
      Caption         =   "Juice Probability:"
      Height          =   255
      Left            =   240
      TabIndex        =   37
      Top             =   6840
      Visible         =   0   'False
      Width           =   1305
   End
   Begin VB.Label Label20 
      Caption         =   "msec."
      ForeColor       =   &H000000FF&
      Height          =   255
      Left            =   6600
      TabIndex        =   35
      Top             =   7350
      Visible         =   0   'False
      Width           =   405
   End
   Begin VB.Label Label19 
      Caption         =   "sec."
      Height          =   270
      Left            =   1920
      TabIndex        =   33
      Top             =   3600
      Width           =   375
   End
   Begin VB.Label Label18 
      Caption         =   "msec"
      Height          =   255
      Left            =   2340
      TabIndex        =   32
      Top             =   5760
      Width           =   405
   End
   Begin VB.Label Label17 
      Caption         =   "msec."
      Height          =   255
      Left            =   2400
      TabIndex        =   30
      Top             =   3960
      Visible         =   0   'False
      Width           =   495
   End
   Begin VB.Label Label2 
      Caption         =   "Subject:"
      Height          =   255
      Left            =   120
      TabIndex        =   29
      Top             =   120
      Width           =   615
   End
   Begin VB.Label Label3 
      Caption         =   "Session ID:"
      Height          =   255
      Left            =   120
      TabIndex        =   28
      Top             =   600
      Width           =   855
   End
   Begin VB.Label Label4 
      Caption         =   "Num of Trials (if this isn't 96, bad things happen)"
      Height          =   495
      Left            =   120
      TabIndex        =   27
      Top             =   1440
      Width           =   2175
   End
   Begin VB.Label Label5 
      Alignment       =   1  'Right Justify
      Caption         =   "Probe Aborts Enabled:"
      Height          =   225
      Left            =   2520
      TabIndex        =   26
      Top             =   7170
      Visible         =   0   'False
      Width           =   1605
   End
   Begin VB.Label Label6 
      Caption         =   "Correction Procedure:"
      Height          =   255
      Left            =   120
      TabIndex        =   25
      Top             =   4920
      Width           =   1815
   End
   Begin VB.Label Label7 
      Alignment       =   1  'Right Justify
      Caption         =   "Abort LTO:"
      Height          =   255
      Left            =   0
      TabIndex        =   24
      Top             =   7680
      Visible         =   0   'False
      Width           =   795
   End
   Begin VB.Label Label8 
      Caption         =   "Time Out:"
      Height          =   255
      Left            =   120
      TabIndex        =   23
      Top             =   2760
      Width           =   1095
   End
   Begin VB.Label Label9 
      Caption         =   "List Length:"
      Height          =   255
      Left            =   7320
      TabIndex        =   22
      Top             =   6720
      Visible         =   0   'False
      Width           =   825
   End
   Begin VB.Label Label10 
      Caption         =   "View Time:"
      Height          =   255
      Left            =   5220
      TabIndex        =   21
      Top             =   7320
      Visible         =   0   'False
      Width           =   885
   End
   Begin VB.Label Label11 
      Caption         =   "ITI (min, max). Make equal for single ITI"
      Height          =   615
      Left            =   120
      TabIndex        =   20
      Top             =   3120
      Width           =   2535
   End
   Begin VB.Label Label12 
      Caption         =   "ISI:"
      Height          =   255
      Left            =   120
      TabIndex        =   19
      Top             =   6630
      Visible         =   0   'False
      Width           =   225
   End
   Begin VB.Label Label13 
      Caption         =   "Delay:"
      Height          =   255
      Left            =   120
      TabIndex        =   18
      Top             =   3960
      Visible         =   0   'False
      Width           =   1215
   End
   Begin VB.Label Label15 
      Caption         =   "sec."
      Height          =   255
      Left            =   2220
      TabIndex        =   16
      Top             =   7680
      Visible         =   0   'False
      Width           =   315
   End
   Begin VB.Label Label16 
      Caption         =   "sec."
      Height          =   255
      Left            =   2400
      TabIndex        =   15
      Top             =   2760
      Width           =   375
   End
End
Attribute VB_Name = "frmDefaults"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False









Option Explicit


'the combo box doesn't get "selected"
'when read from a DB apparently
Private combo_changed As Boolean

Private Sub cmdCancel_Click()

    Unload Me
    
End Sub


Private Sub cmdNormalMode_Click()

    lblOpsMode.Caption = "Normal"
    
End Sub

Private Sub cmdOK_Click()
    
    Dim bError As Boolean
    
    bError = False
    
    Dim i, s As Integer
    s = 0
    For i = 1 To 6
       s = s + VT(i)
    Next
    
    For i = 12 To 16
       s = s + VT(i)
    Next
    
    If s <> intNumTrials Then
        bError = True
        MsgBox "Sum of display sizes MUST equal Trial Count: " & s & " <> " & intNumTrials, vbExclamation, "Abort Corrections"
    End If
    
    If Int(Combo2.Text) < 8 Then
        s = 0
        For i = 4 To 6
           s = s + VT(i)
        Next
        
        If (s > 0) Then
            MsgBox "Set size of " & Combo2.Text & " is incompatible with desired display size counts", vbCritical
            bError = True
        End If
    End If
    
    
    If chkAbortCorrections = 1 And intAbortLTO = 0 Then
        bError = True
        MsgBox "When Abort Corrections are enabled, LTO can not be 0", vbExclamation, "Abort Corrections"
    End If
   
    If chkIncorrectCorrections = 1 And intIncorrectLTO = 0 Then
        bError = True
        MsgBox "When Incorrect Corrections are enabled, LTO can not be 0", vbExclamation, "Incorrect Corrections"
    End If
 
    If bError = False Then
        
        gsImagePath = dirImages.Path
        datDefaults.Recordset.Edit
        datDefaults.Recordset("ImagePath") = gsImagePath
        datDefaults.Recordset.Update
        
        frmMain.mnuSetup.Enabled = True
        frmMain.mnuRun.Enabled = True
        frmMain.mnuResults.Enabled = True
        Unload Me
        
    End If
    
End Sub



Private Sub cmdSoundFile_Click()
    
    On Error GoTo Error_cmdSoundFile_Click
    
    diaSound.Filter = "Sound File (*.wav)|*.wav|All Files (*.*)|*.*"
    diaSound.DialogTitle = "Open Sound File"
    diaSound.DefaultExt = "wav"
    diaSound.InitDir = App.Path
    diaSound.CancelError = True
    diaSound.Flags = cdlOFNFileMustExist + cdlOFNHideReadOnly
    diaSound.ShowOpen
    
    lblSound.Caption = diaSound.FileName
    
    If lblSound.Caption <> "" Then
        cmdTestSound.Enabled = True
    Else
        cmdTestSound.Enabled = False
    End If
    
    Exit Sub
    
Error_cmdSoundFile_Click:

    If Err = 76 Then  'file error
        Resume Next
    Else
        Exit Sub
    End If
    

End Sub

Private Sub cmdTestMode_Click()

    lblOpsMode.Caption = "Test"
    
End Sub


Private Sub cmdTestSound_Click()

    On Error GoTo Error_cmdTestSound_Click
    
    If lblSound.Caption <> "" Then
        waveTest.WaveIndex = 0
        waveTest.FilenameRead = lblSound.Caption
        waveTest.Play = 0
    End If
    
    Exit Sub
    
Error_cmdTestSound_Click:

    MsgBox Err & ":  " & Error
    Exit Sub
    
End Sub



Private Sub Combo2_Change()
    combo_changed = True
End Sub


Private Sub Form_Activate()

    Refresh
    
    datDefaults.DatabaseName = gsProfile
    datDefaults.RecordSource = "SessionInfo"
    datDefaults.Refresh
    
    datView.DatabaseName = gsProfile
    datView.RecordSource = "ViewConstants"
    datView.Refresh

    optPhase(datDefaults.Recordset("Phase") - 1).Value = True
    
    lblDate.Caption = Format(Now, "yy/mm/dd hh:mm:ss")
    txtSessionFile.Text = gsProfile
    
    On Error GoTo Error_ImagePath
    dirImages.Path = datDefaults.Recordset("ImagePath")
    lblResultsFile.Caption = gsResultsFile

    If lblSound.Caption <> "" Then
        cmdTestSound.Enabled = True
    Else
        cmdTestSound.Enabled = False
    End If

    Exit Sub
    
Error_ImagePath:
    dirImages.Path = "c:\"
    Resume Next
    
End Sub


Private Sub Form_Load()
    
    gsubScreenDimensions Me, 0, 0
        combo_changed = False
End Sub


Private Sub optPhase_Click(Index As Integer)

    Select Case Index
    
        Case 0:
            lblPhase.Caption = 1
    
        Case 1:
            lblPhase.Caption = 2
    
        Case 2:
            lblPhase.Caption = 3
    
    End Select
    
End Sub


