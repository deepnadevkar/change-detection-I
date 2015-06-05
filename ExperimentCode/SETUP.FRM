VERSION 5.00
Object = "{02B97C13-ED1E-11CD-A08B-00AA00575482}#1.0#0"; "mhiinp32.ocx"
Object = "{3AE5AE83-A6DA-101B-9313-00AA00575482}#1.0#0"; "Mhfram32.ocx"
Object = "{8D650141-6025-11D1-BC40-0000C042AEC0}#3.0#0"; "ssdw3b32.ocx"
Begin VB.Form frmSetup 
   BackColor       =   &H00404040&
   BorderStyle     =   1  'Fixed Single
   Caption         =   "Setup"
   ClientHeight    =   8025
   ClientLeft      =   810
   ClientTop       =   1995
   ClientWidth     =   11475
   ControlBox      =   0   'False
   LinkTopic       =   "Form1"
   MaxButton       =   0   'False
   MinButton       =   0   'False
   PaletteMode     =   1  'UseZOrder
   ScaleHeight     =   8025
   ScaleWidth      =   11475
   Begin VB.Data datProbeInfo 
      Caption         =   "ProbeInfo"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   345
      Left            =   1560
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   4680
      Visible         =   0   'False
      Width           =   2025
   End
   Begin VB.Data datView 
      Caption         =   "View"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   5820
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   7440
      Visible         =   0   'False
      Width           =   2145
   End
   Begin VB.TextBox ErrorBox 
      Alignment       =   2  'Center
      Appearance      =   0  'Flat
      BackColor       =   &H000000FF&
      BeginProperty Font 
         Name            =   "MS Sans Serif"
         Size            =   18
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      ForeColor       =   &H00FFFFFF&
      Height          =   795
      Left            =   4680
      TabIndex        =   83
      Text            =   "Text1"
      Top             =   6840
      Visible         =   0   'False
      Width           =   4215
   End
   Begin Mh3dfrmLibCtl.Mh3dFrame frameStatus 
      Height          =   645
      Left            =   3240
      TabIndex        =   80
      Top             =   2400
      Visible         =   0   'False
      Width           =   5325
      _Version        =   65536
      _ExtentX        =   9393
      _ExtentY        =   1138
      _StockProps     =   79
      Caption         =   "xxx"
      BeginProperty Font {0BE35203-8F91-11CE-9DE3-00AA004BB851} 
         Name            =   "Arial"
         Size            =   14.25
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Alignment       =   2
      AutoSize        =   0   'False
      BevelSize       =   3
      BevelStyle      =   1
      BorderColor     =   -2147483642
      BorderStyle     =   1
      FillColor       =   12632064
      FontStyle       =   0
      FontTransparent =   0   'False
      LightColor      =   -2147483643
      ShadowColor     =   -2147483632
      TextColor       =   -2147483640
      WallPaper       =   0
      NoPrefix        =   0   'False
      FormatString    =   ""
      Caption         =   "xxx"
      Picture         =   "SETUP.frx":0000
   End
   Begin Mh3dfrmLibCtl.Mh3dFrame frameMessage 
      Height          =   375
      Left            =   3750
      TabIndex        =   79
      Top             =   1830
      Visible         =   0   'False
      Width           =   2205
      _Version        =   65536
      _ExtentX        =   3889
      _ExtentY        =   661
      _StockProps     =   79
      Caption         =   "Greater than List Length."
      BeginProperty Font {0BE35203-8F91-11CE-9DE3-00AA004BB851} 
         Name            =   "MS Sans Serif"
         Size            =   8.25
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Alignment       =   2
      AutoSize        =   0   'False
      BevelSize       =   2
      BevelStyle      =   1
      BorderColor     =   -2147483642
      BorderStyle     =   1
      FillColor       =   192
      FontStyle       =   0
      FontTransparent =   0   'False
      LightColor      =   -2147483643
      ShadowColor     =   -2147483632
      TextColor       =   16777215
      WallPaper       =   0
      NoPrefix        =   0   'False
      FormatString    =   ""
      Caption         =   "Greater than List Length."
      Picture         =   "SETUP.frx":001C
   End
   Begin VB.CommandButton cmdNextReinforceNoapparatus 
      Caption         =   "Next Reinf, No App"
      Height          =   345
      Left            =   7140
      TabIndex        =   75
      Top             =   1470
      Width           =   1635
   End
   Begin VB.CommandButton cmdNextSame0Position 
      Caption         =   "Next Same, 0 Pos"
      Height          =   345
      Left            =   5460
      TabIndex        =   74
      Top             =   1470
      Width           =   1665
   End
   Begin VB.Frame Frame9 
      Height          =   945
      Left            =   9240
      TabIndex        =   16
      Top             =   0
      Width           =   2235
      Begin VB.CommandButton cmdJuice 
         Caption         =   "Juice"
         Height          =   255
         Left            =   1200
         TabIndex        =   35
         Top             =   330
         Width           =   765
      End
      Begin VB.CommandButton cmdPellets 
         Caption         =   "Pellets"
         Height          =   255
         Left            =   1200
         TabIndex        =   34
         Top             =   630
         Width           =   765
      End
      Begin VB.Label lblApparatus 
         Alignment       =   2  'Center
         BackColor       =   &H00FFFFFF&
         BorderStyle     =   1  'Fixed Single
         Caption         =   "Pellets"
         DataField       =   "Apparatus"
         DataSource      =   "datTrialInfo"
         Height          =   255
         Left            =   150
         TabIndex        =   78
         Top             =   600
         Width           =   765
      End
      Begin VB.Label Label10 
         Caption         =   "Apparatus:"
         Height          =   255
         Left            =   150
         TabIndex        =   36
         Top             =   360
         Width           =   825
      End
   End
   Begin VB.Frame Frame1 
      Height          =   945
      Left            =   6420
      TabIndex        =   31
      Top             =   0
      Width           =   2835
      Begin VB.CheckBox chkNoTouch 
         DataField       =   "NoTouchTrial"
         DataSource      =   "datTrialInfo"
         Height          =   495
         Left            =   960
         TabIndex        =   84
         Top             =   480
         Width           =   255
      End
      Begin VB.CheckBox chkTransfer 
         DataField       =   "Transfer"
         DataSource      =   "datTrialInfo"
         Height          =   495
         Left            =   2280
         TabIndex        =   81
         Top             =   120
         Width           =   255
      End
      Begin VB.CheckBox chkReinforce 
         DataField       =   "Reinforce"
         DataSource      =   "datTrialInfo"
         Height          =   495
         Left            =   960
         TabIndex        =   33
         Top             =   120
         Width           =   255
      End
      Begin VB.Label Label18 
         Caption         =   "No Peck:"
         Height          =   255
         Left            =   120
         TabIndex        =   85
         Top             =   600
         Width           =   855
      End
      Begin VB.Label Label17 
         Caption         =   "Transfer:"
         Height          =   255
         Left            =   1560
         TabIndex        =   82
         Top             =   240
         Width           =   975
      End
      Begin VB.Label Label8 
         Caption         =   "Reinforce:"
         Height          =   255
         Left            =   120
         TabIndex        =   32
         Top             =   240
         Width           =   855
      End
   End
   Begin VB.Frame Frame8 
      Height          =   945
      Left            =   5490
      TabIndex        =   15
      Top             =   0
      Width           =   945
      Begin MhinintLibCtl.MhIntInput intITI 
         DataField       =   "ITI"
         DataSource      =   "datTrialInfo"
         Height          =   285
         Left            =   60
         TabIndex        =   29
         Top             =   600
         Width           =   405
         _Version        =   65536
         _ExtentX        =   714
         _ExtentY        =   503
         _StockProps     =   77
         BackColor       =   -2147483633
         BevelSize       =   1
         BorderStyle     =   0
         FillColor       =   16777215
         Min             =   -2147483647
         AutoHScroll     =   -1  'True
         CaretColor      =   -2147483642
      End
      Begin VB.Label Label5 
         Caption         =   "ITI:"
         Height          =   225
         Left            =   300
         TabIndex        =   30
         Top             =   360
         Width           =   315
      End
      Begin VB.Label Label11 
         Caption         =   "sec."
         Height          =   255
         Left            =   510
         TabIndex        =   28
         Top             =   630
         Width           =   345
      End
   End
   Begin VB.Frame Frame7 
      Height          =   945
      Left            =   4680
      TabIndex        =   14
      Top             =   0
      Width           =   825
      Begin MhinintLibCtl.MhIntInput intListLength 
         DataField       =   "ListLength"
         DataSource      =   "datTrialInfo"
         Height          =   285
         Left            =   210
         TabIndex        =   26
         Top             =   600
         Width           =   405
         _Version        =   65536
         _ExtentX        =   714
         _ExtentY        =   503
         _StockProps     =   77
         BackColor       =   -2147483633
         BevelSize       =   1
         BorderStyle     =   0
         FillColor       =   12632256
         Min             =   -2147483647
         AutoHScroll     =   -1  'True
         CaretColor      =   -2147483642
      End
      Begin VB.Label Label4 
         Caption         =   "List Length:"
         Height          =   405
         Left            =   150
         TabIndex        =   27
         Top             =   180
         Width           =   615
      End
   End
   Begin VB.Frame Frame6 
      Height          =   945
      Left            =   3810
      TabIndex        =   13
      Top             =   0
      Width           =   885
      Begin MhinintLibCtl.MhIntInput intCorrectPosition 
         DataField       =   "CorrectPosition"
         DataSource      =   "datTrialInfo"
         Height          =   285
         Left            =   210
         TabIndex        =   24
         Top             =   600
         Width           =   405
         _Version        =   65536
         _ExtentX        =   714
         _ExtentY        =   503
         _StockProps     =   77
         BackColor       =   16777215
         BevelSize       =   1
         BorderStyle     =   0
         FillColor       =   16777215
         Min             =   0
         AutoHScroll     =   -1  'True
         CaretColor      =   -2147483642
      End
      Begin VB.Label Label3 
         Caption         =   "Correct Position:"
         Height          =   435
         Left            =   150
         TabIndex        =   25
         Top             =   180
         Width           =   675
      End
   End
   Begin VB.Frame Frame5 
      Height          =   945
      Left            =   1950
      TabIndex        =   12
      Top             =   0
      Width           =   1875
      Begin VB.CommandButton cmdSame 
         Caption         =   "Left"
         Height          =   255
         Left            =   1020
         TabIndex        =   22
         Top             =   330
         Width           =   765
      End
      Begin VB.CommandButton cmdDiff 
         Caption         =   "Right"
         Height          =   255
         Left            =   1020
         TabIndex        =   21
         Top             =   630
         Width           =   765
      End
      Begin VB.Label lblType 
         Alignment       =   2  'Center
         BackColor       =   &H00FFFFFF&
         BorderStyle     =   1  'Fixed Single
         Caption         =   "Left"
         DataField       =   "TrialType"
         DataSource      =   "datTrialInfo"
         Height          =   255
         Left            =   120
         TabIndex        =   77
         Top             =   630
         Width           =   705
      End
      Begin VB.Label Label2 
         Caption         =   "Trial Type:"
         Height          =   255
         Left            =   90
         TabIndex        =   23
         Top             =   360
         Width           =   825
      End
   End
   Begin VB.Frame Frame4 
      Height          =   945
      Left            =   0
      TabIndex        =   11
      Top             =   0
      Width           =   1965
      Begin VB.Label Label1 
         Caption         =   "Trial#:"
         Height          =   255
         Left            =   90
         TabIndex        =   20
         Top             =   420
         Width           =   555
      End
      Begin VB.Label lblTrialNum 
         Alignment       =   2  'Center
         Caption         =   "000"
         DataField       =   "TrialNum"
         DataSource      =   "datTrialInfo"
         BeginProperty Font 
            Name            =   "MS Sans Serif"
            Size            =   13.5
            Charset         =   0
            Weight          =   400
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   375
         Left            =   630
         TabIndex        =   19
         Top             =   330
         Width           =   465
      End
      Begin VB.Label Label9 
         Caption         =   "of"
         Height          =   225
         Left            =   1170
         TabIndex        =   18
         Top             =   420
         Width           =   225
      End
      Begin VB.Label lblNumTrials 
         Alignment       =   2  'Center
         Caption         =   "000"
         DataField       =   "NumOfTrials"
         DataSource      =   "datSessionInfo"
         BeginProperty Font 
            Name            =   "MS Sans Serif"
            Size            =   13.5
            Charset         =   0
            Weight          =   400
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   375
         Left            =   1380
         TabIndex        =   17
         Top             =   330
         Width           =   465
      End
   End
   Begin VB.Data datStats 
      Caption         =   "Stats"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   2730
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   7350
      Visible         =   0   'False
      Width           =   1755
   End
   Begin VB.Data datSames 
      Caption         =   "Sames"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   2730
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   7050
      Visible         =   0   'False
      Width           =   1815
   End
   Begin VB.Data datRecordCounter 
      Caption         =   "RecordCounter"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   2730
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   6750
      Visible         =   0   'False
      Width           =   2415
   End
   Begin VB.Data datSessionInfo 
      Caption         =   "SessionInfo"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   2940
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   8010
      Visible         =   0   'False
      Width           =   2145
   End
   Begin VB.Data datTrialInfo 
      BackColor       =   &H0000C000&
      Caption         =   "Trial Selection"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      BeginProperty Font 
         Name            =   "MS Sans Serif"
         Size            =   13.5
         Charset         =   0
         Weight          =   700
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Height          =   495
      Left            =   5460
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   960
      Width           =   3315
   End
   Begin VB.Data datProbe 
      Caption         =   "Probe"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   5820
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   7140
      Visible         =   0   'False
      Width           =   2025
   End
   Begin VB.Data datStimuli 
      Caption         =   "Stimuli"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   345
      Left            =   5820
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   6810
      Visible         =   0   'False
      Width           =   2445
   End
   Begin VB.Frame Frame2 
      Height          =   5745
      Left            =   0
      TabIndex        =   2
      Top             =   840
      Width           =   5445
      Begin VB.FileListBox fileImages 
         Height          =   285
         Left            =   480
         Pattern         =   "*.tga"
         TabIndex        =   3
         Top             =   4890
         Visible         =   0   'False
         Width           =   2085
      End
      Begin SSDataWidgets_B.SSDBGrid gridTrialInfo 
         Bindings        =   "SETUP.frx":0038
         Height          =   4935
         Left            =   120
         TabIndex        =   76
         Top             =   600
         Width           =   5145
         _Version        =   196617
         RecordSelectors =   0   'False
         AllowAddNew     =   -1  'True
         AllowDelete     =   -1  'True
         AllowRowSizing  =   0   'False
         AllowGroupSizing=   0   'False
         AllowColumnSizing=   0   'False
         AllowGroupMoving=   0   'False
         AllowGroupSwapping=   0   'False
         AllowColumnSwapping=   0
         AllowGroupShrinking=   0   'False
         AllowColumnShrinking=   0   'False
         AllowDragDrop   =   0   'False
         SelectTypeRow   =   1
         ForeColorEven   =   0
         BackColorOdd    =   16777215
         RowHeight       =   423
         Columns.Count   =   8
         Columns(0).Width=   582
         Columns(0).Caption=   "#"
         Columns(0).Name =   "#"
         Columns(0).Alignment=   1
         Columns(0).CaptionAlignment=   1
         Columns(0).AllowSizing=   0   'False
         Columns(0).DataField=   "TrialNum"
         Columns(0).DataType=   3
         Columns(0).FieldLen=   256
         Columns(0).Locked=   -1  'True
         Columns(1).Width=   873
         Columns(1).Caption=   "Type"
         Columns(1).Name =   "Type"
         Columns(1).CaptionAlignment=   0
         Columns(1).DataField=   "TrialType"
         Columns(1).DataType=   8
         Columns(1).FieldLen=   256
         Columns(2).Width=   1984
         Columns(2).Caption=   "Configuration"
         Columns(2).Name =   "configuration"
         Columns(2).DataField=   "Configuration"
         Columns(2).FieldLen=   256
         Columns(3).Width=   900
         Columns(3).Caption=   "Reinf"
         Columns(3).Name =   "Reinf"
         Columns(3).Alignment=   1
         Columns(3).CaptionAlignment=   1
         Columns(3).DataField=   "Reinforce"
         Columns(3).DataType=   11
         Columns(3).FieldLen=   256
         Columns(3).Style=   2
         Columns(4).Width=   979
         Columns(4).Caption=   "App"
         Columns(4).Name =   "App"
         Columns(4).Alignment=   2
         Columns(4).CaptionAlignment=   2
         Columns(4).DataField=   "Apparatus"
         Columns(4).DataType=   8
         Columns(4).FieldLen=   256
         Columns(5).Width=   1429
         Columns(5).Caption=   "NoTouch"
         Columns(5).Name =   "NoTouch"
         Columns(5).DataField=   "NoTouchTrial"
         Columns(5).FieldLen=   256
         Columns(5).Style=   2
         Columns(6).Width=   582
         Columns(6).Caption=   "ITI"
         Columns(6).Name =   "ITI"
         Columns(6).Alignment=   1
         Columns(6).CaptionAlignment=   1
         Columns(6).DataField=   "ITI"
         Columns(6).DataType=   3
         Columns(6).FieldLen=   256
         Columns(7).Width=   820
         Columns(7).Caption=   "Tran"
         Columns(7).Name =   "Tran"
         Columns(7).CaptionAlignment=   0
         Columns(7).DataField=   "Transfer"
         Columns(7).DataType=   11
         Columns(7).FieldLen=   256
         Columns(7).Style=   2
         UseDefaults     =   0   'False
         _ExtentX        =   9075
         _ExtentY        =   8705
         _StockProps     =   79
         Caption         =   "Trial Info"
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
      Begin VB.Label lblImagePath 
         BorderStyle     =   1  'Fixed Single
         DataField       =   "ImagePath"
         DataSource      =   "datSessionInfo"
         Height          =   285
         Left            =   960
         TabIndex        =   10
         Top             =   240
         Width           =   3315
      End
   End
   Begin VB.Frame Frame11 
      Height          =   855
      Left            =   8760
      TabIndex        =   49
      Top             =   960
      Width           =   2715
      Begin VB.Label Label31 
         Caption         =   "Reinforce"
         Height          =   225
         Left            =   60
         TabIndex        =   58
         Top             =   510
         Width           =   705
      End
      Begin VB.Label lblReinforceCount 
         Alignment       =   1  'Right Justify
         Caption         =   "00"
         BeginProperty Font 
            Name            =   "MS Sans Serif"
            Size            =   8.25
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   225
         Left            =   780
         TabIndex        =   57
         Top             =   510
         Width           =   315
      End
      Begin VB.Label Label29 
         Caption         =   "Count"
         Height          =   255
         Left            =   750
         TabIndex        =   56
         Top             =   210
         Width           =   495
      End
      Begin VB.Label Label28 
         Caption         =   "Actual"
         Height          =   255
         Left            =   1350
         TabIndex        =   55
         Top             =   210
         Width           =   495
      End
      Begin VB.Label Label27 
         Caption         =   "Default"
         Height          =   255
         Left            =   2010
         TabIndex        =   54
         Top             =   210
         Width           =   555
      End
      Begin VB.Label Label26 
         Caption         =   "%"
         Height          =   225
         Left            =   1680
         TabIndex        =   53
         Top             =   510
         Width           =   195
      End
      Begin VB.Label lblReinforceActual 
         Alignment       =   1  'Right Justify
         Caption         =   "00"
         BeginProperty Font 
            Name            =   "MS Sans Serif"
            Size            =   8.25
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   225
         Left            =   1320
         TabIndex        =   52
         Top             =   510
         Width           =   315
      End
      Begin VB.Label lblReinforceDefault 
         Alignment       =   1  'Right Justify
         Caption         =   "00"
         DataField       =   "ReinfProbability"
         DataSource      =   "datSessionInfo"
         BeginProperty Font 
            Name            =   "MS Sans Serif"
            Size            =   8.25
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   225
         Left            =   2040
         TabIndex        =   51
         Top             =   510
         Width           =   315
      End
      Begin VB.Label Label23 
         Caption         =   "%"
         Height          =   225
         Left            =   2400
         TabIndex        =   50
         Top             =   510
         Width           =   195
      End
   End
   Begin VB.Frame Frame10 
      Height          =   1125
      Left            =   8760
      TabIndex        =   48
      Top             =   1710
      Width           =   2715
      Begin VB.Label Label46 
         Caption         =   "Juice"
         Height          =   225
         Left            =   240
         TabIndex        =   73
         Top             =   510
         Visible         =   0   'False
         Width           =   555
      End
      Begin VB.Label Label45 
         Caption         =   "Pellets"
         Height          =   225
         Left            =   210
         TabIndex        =   72
         Top             =   810
         Width           =   555
      End
      Begin VB.Label lblJuiceCount 
         Alignment       =   1  'Right Justify
         Caption         =   "00"
         BeginProperty Font 
            Name            =   "MS Sans Serif"
            Size            =   8.25
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   225
         Left            =   810
         TabIndex        =   71
         Top             =   510
         Visible         =   0   'False
         Width           =   315
      End
      Begin VB.Label lblPelletsCount 
         Alignment       =   1  'Right Justify
         Caption         =   "00"
         BeginProperty Font 
            Name            =   "MS Sans Serif"
            Size            =   8.25
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   225
         Left            =   780
         TabIndex        =   70
         Top             =   810
         Width           =   315
      End
      Begin VB.Label Label42 
         Caption         =   "Count"
         Height          =   255
         Left            =   750
         TabIndex        =   69
         Top             =   210
         Width           =   495
      End
      Begin VB.Label Label41 
         Caption         =   "Actual"
         Height          =   255
         Left            =   1350
         TabIndex        =   68
         Top             =   210
         Width           =   495
      End
      Begin VB.Label Label40 
         Caption         =   "Default"
         Height          =   255
         Left            =   2010
         TabIndex        =   67
         Top             =   210
         Width           =   555
      End
      Begin VB.Label Label39 
         Caption         =   "%"
         Height          =   225
         Left            =   1710
         TabIndex        =   66
         Top             =   510
         Visible         =   0   'False
         Width           =   195
      End
      Begin VB.Label lblPelletsActual 
         Alignment       =   1  'Right Justify
         Caption         =   "00"
         BeginProperty Font 
            Name            =   "MS Sans Serif"
            Size            =   8.25
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   225
         Left            =   1320
         TabIndex        =   65
         Top             =   810
         Width           =   315
      End
      Begin VB.Label lblJuiceActual 
         Alignment       =   1  'Right Justify
         Caption         =   "00"
         BeginProperty Font 
            Name            =   "MS Sans Serif"
            Size            =   8.25
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   225
         Left            =   1350
         TabIndex        =   64
         Top             =   510
         Visible         =   0   'False
         Width           =   315
      End
      Begin VB.Label lblPelletsDefault 
         Alignment       =   1  'Right Justify
         Caption         =   "00"
         DataField       =   "PelletProbability"
         DataSource      =   "datSessionInfo"
         BeginProperty Font 
            Name            =   "MS Sans Serif"
            Size            =   8.25
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   225
         Left            =   2040
         TabIndex        =   63
         Top             =   810
         Width           =   315
      End
      Begin VB.Label lblJuiceDefault 
         Alignment       =   1  'Right Justify
         Caption         =   "00"
         DataField       =   "JuiceProbability"
         DataSource      =   "datSessionInfo"
         BeginProperty Font 
            Name            =   "MS Sans Serif"
            Size            =   8.25
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   225
         Left            =   2070
         TabIndex        =   62
         Top             =   510
         Visible         =   0   'False
         Width           =   315
      End
      Begin VB.Label Label34 
         Caption         =   "%"
         Height          =   225
         Left            =   1680
         TabIndex        =   61
         Top             =   810
         Width           =   195
      End
      Begin VB.Label Label33 
         Caption         =   "%"
         Height          =   225
         Left            =   2430
         TabIndex        =   60
         Top             =   510
         Visible         =   0   'False
         Width           =   195
      End
      Begin VB.Label Label32 
         Caption         =   "%"
         Height          =   225
         Left            =   2400
         TabIndex        =   59
         Top             =   810
         Width           =   195
      End
   End
   Begin VB.Frame Frame3 
      Height          =   4215
      Left            =   8760
      TabIndex        =   4
      Top             =   2730
      Width           =   2715
      Begin SSDataWidgets_B.SSDBGrid gridSames 
         Bindings        =   "SETUP.frx":0053
         Height          =   2925
         Left            =   180
         TabIndex        =   7
         Top             =   1170
         Width           =   2445
         _Version        =   196617
         AllowAddNew     =   -1  'True
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
         SelectTypeCol   =   0
         SelectTypeRow   =   0
         ForeColorEven   =   0
         BackColorOdd    =   16777215
         RowHeight       =   423
         Columns(0).Width=   3200
         UseDefaults     =   0   'False
         _ExtentX        =   4313
         _ExtentY        =   5159
         _StockProps     =   79
         Caption         =   "Statistics"
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
      Begin VB.Label Label22 
         Caption         =   "%"
         Height          =   225
         Left            =   2400
         TabIndex        =   47
         Top             =   780
         Width           =   195
      End
      Begin VB.Label Label21 
         Caption         =   "%"
         Height          =   225
         Left            =   2400
         TabIndex        =   46
         Top             =   480
         Width           =   195
      End
      Begin VB.Label Label20 
         Caption         =   "%"
         Height          =   225
         Left            =   1680
         TabIndex        =   45
         Top             =   780
         Width           =   195
      End
      Begin VB.Label lblDiffsDefault 
         Alignment       =   1  'Right Justify
         Caption         =   "00"
         DataField       =   "DiffProbability"
         DataSource      =   "datSessionInfo"
         BeginProperty Font 
            Name            =   "MS Sans Serif"
            Size            =   8.25
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   225
         Left            =   2040
         TabIndex        =   44
         Top             =   480
         Width           =   315
      End
      Begin VB.Label lblSamesDefault 
         Alignment       =   1  'Right Justify
         Caption         =   "00"
         DataField       =   "SameProbability"
         DataSource      =   "datSessionInfo"
         BeginProperty Font 
            Name            =   "MS Sans Serif"
            Size            =   8.25
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   225
         Left            =   2040
         TabIndex        =   43
         Top             =   780
         Width           =   315
      End
      Begin VB.Label lblDiffsActual 
         Alignment       =   1  'Right Justify
         Caption         =   "00"
         BeginProperty Font 
            Name            =   "MS Sans Serif"
            Size            =   8.25
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   225
         Left            =   1320
         TabIndex        =   42
         Top             =   480
         Width           =   315
      End
      Begin VB.Label lblSamesActual 
         Alignment       =   1  'Right Justify
         Caption         =   "00"
         BeginProperty Font 
            Name            =   "MS Sans Serif"
            Size            =   8.25
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   225
         Left            =   1320
         TabIndex        =   41
         Top             =   780
         Width           =   315
      End
      Begin VB.Label Label15 
         Caption         =   "%"
         Height          =   225
         Left            =   1680
         TabIndex        =   40
         Top             =   480
         Width           =   195
      End
      Begin VB.Label Label14 
         Caption         =   "Default"
         Height          =   255
         Left            =   2040
         TabIndex        =   39
         Top             =   180
         Width           =   555
      End
      Begin VB.Label Label13 
         Caption         =   "Actual"
         Height          =   255
         Left            =   1380
         TabIndex        =   38
         Top             =   180
         Width           =   495
      End
      Begin VB.Label Label12 
         Caption         =   "Count"
         Height          =   255
         Left            =   780
         TabIndex        =   37
         Top             =   180
         Width           =   495
      End
      Begin VB.Label lblSamesCount 
         Alignment       =   1  'Right Justify
         Caption         =   "00"
         BeginProperty Font 
            Name            =   "MS Sans Serif"
            Size            =   8.25
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   225
         Left            =   780
         TabIndex        =   9
         Top             =   780
         Width           =   315
      End
      Begin VB.Label lblDiffsCount 
         Alignment       =   1  'Right Justify
         Caption         =   "00"
         BeginProperty Font 
            Name            =   "MS Sans Serif"
            Size            =   8.25
            Charset         =   0
            Weight          =   700
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         Height          =   225
         Left            =   780
         TabIndex        =   8
         Top             =   480
         Width           =   315
      End
      Begin VB.Label Label7 
         Caption         =   "Lefts"
         Height          =   225
         Left            =   210
         TabIndex        =   6
         Top             =   780
         Width           =   555
      End
      Begin VB.Label Label6 
         Caption         =   "Rights"
         Height          =   225
         Left            =   210
         TabIndex        =   5
         Top             =   480
         Width           =   555
      End
   End
   Begin SSDataWidgets_B.SSDBGrid gridProbe 
      Bindings        =   "SETUP.frx":006A
      Height          =   1095
      Left            =   5400
      TabIndex        =   1
      Top             =   3720
      Width           =   3405
      _Version        =   196617
      SelectTypeRow   =   1
      ForeColorEven   =   0
      BackColorOdd    =   16777215
      RowHeight       =   423
      Columns.Count   =   3
      Columns(0).Width=   2381
      Columns(0).Caption=   "Filename"
      Columns(0).CaptionAlignment=   0
      Columns(0).DataField=   "Filename"
      Columns(0).DataType=   8
      Columns(0).FieldLen=   256
      Columns(1).Width=   1164
      Columns(1).Caption=   "Time"
      Columns(1).Alignment=   1
      Columns(1).CaptionAlignment=   1
      Columns(1).DataField=   "ViewTime"
      Columns(1).DataType=   3
      Columns(1).FieldLen=   256
      Columns(2).Width=   1191
      Columns(2).Caption=   "Delay"
      Columns(2).Alignment=   1
      Columns(2).CaptionAlignment=   1
      Columns(2).DataField=   "Delay"
      Columns(2).DataType=   3
      Columns(2).FieldLen=   256
      UseDefaults     =   0   'False
      _ExtentX        =   6006
      _ExtentY        =   1931
      _StockProps     =   79
      Caption         =   "Item 2"
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
   Begin SSDataWidgets_B.SSDBGrid gridStimuli 
      Bindings        =   "SETUP.frx":0081
      Height          =   1395
      Left            =   5400
      TabIndex        =   0
      Top             =   1860
      Width           =   3405
      _Version        =   196617
      SelectTypeCol   =   0
      SelectTypeRow   =   1
      ForeColorEven   =   0
      BackColorOdd    =   16777215
      RowHeight       =   423
      Columns.Count   =   4
      Columns(0).Width=   688
      Columns(0).Caption=   "Pos"
      Columns(0).Alignment=   1
      Columns(0).CaptionAlignment=   1
      Columns(0).DataField=   "Position"
      Columns(0).DataType=   3
      Columns(0).FieldLen=   256
      Columns(0).Locked=   -1  'True
      Columns(0).HasBackColor=   -1  'True
      Columns(0).BackColor=   12632256
      Columns(1).Width=   2223
      Columns(1).Caption=   "Filename"
      Columns(1).CaptionAlignment=   0
      Columns(1).DataField=   "Filename"
      Columns(1).DataType=   8
      Columns(1).FieldLen=   256
      Columns(1).HasBackColor=   -1  'True
      Columns(1).BackColor=   16777215
      Columns(2).Width=   1164
      Columns(2).Caption=   "Time"
      Columns(2).Alignment=   1
      Columns(2).CaptionAlignment=   1
      Columns(2).DataField=   "ViewTime"
      Columns(2).DataType=   3
      Columns(2).FieldLen=   256
      Columns(2).HasBackColor=   -1  'True
      Columns(2).BackColor=   16777215
      Columns(3).Width=   741
      Columns(3).Caption=   "ISI"
      Columns(3).Alignment=   1
      Columns(3).CaptionAlignment=   1
      Columns(3).DataField=   "ISI"
      Columns(3).DataType=   3
      Columns(3).FieldLen=   256
      Columns(3).HasBackColor=   -1  'True
      Columns(3).BackColor=   16777215
      UseDefaults     =   0   'False
      _ExtentX        =   6006
      _ExtentY        =   2461
      _StockProps     =   79
      Caption         =   "Item 1"
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
   Begin SSDataWidgets_B.SSDBGrid SSDBGrid1 
      Bindings        =   "SETUP.frx":009A
      Height          =   1215
      Left            =   5400
      TabIndex        =   86
      Top             =   5040
      Width           =   3405
      _Version        =   196617
      SelectTypeRow   =   1
      ForeColorEven   =   0
      BackColorOdd    =   16777215
      RowHeight       =   423
      Columns.Count   =   3
      Columns(0).Width=   2275
      Columns(0).Caption=   "Filename2"
      Columns(0).Name =   "Filename2"
      Columns(0).DataField=   "Filename2"
      Columns(0).FieldLen=   256
      Columns(1).Width=   1667
      Columns(1).Caption=   "Time"
      Columns(1).Name =   "Time"
      Columns(1).Alignment=   1
      Columns(1).CaptionAlignment=   1
      Columns(1).DataField=   "ViewTime"
      Columns(1).DataType=   3
      Columns(1).FieldLen=   256
      Columns(2).Width=   1191
      Columns(2).Caption=   "Delay"
      Columns(2).Name =   "Delay"
      Columns(2).Alignment=   1
      Columns(2).CaptionAlignment=   1
      Columns(2).DataField=   "Delay"
      Columns(2).DataType=   3
      Columns(2).FieldLen=   256
      UseDefaults     =   0   'False
      _ExtentX        =   6006
      _ExtentY        =   2143
      _StockProps     =   79
      Caption         =   "Item 3"
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
   Begin VB.Menu mnuCancelSetup 
      Caption         =   "Cancel &Setup"
   End
   Begin VB.Menu mnuUpdateProbabilities 
      Caption         =   "&Update Probabilities"
   End
End
Attribute VB_Name = "frmSetup"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False
Option Explicit

Dim bLoadingForm As Boolean
Dim sActiveGrid As String
Dim bLoadingRecords As Boolean
Dim bPrinting As Boolean
Dim StimName As String
Dim PickBuffS(240) As Integer
Dim PickBuffD(168) As Integer

Private Sub PickWoReplaceS(Size As Integer, Pick)
      
    Dim Flag As Integer
    Dim CHECK As Integer
    
    Randomize
    Flag = 0
    For CHECK = 1 To Size
        If PickBuffS(CHECK) = 0 Then
            Flag = 1
        End If
    Next CHECK
    
    If Flag = 1 Then
AGAIN:
        Pick = Int(Rnd * Size) + 1
        If PickBuffS(Pick) = 1 Then
            GoTo AGAIN
        End If
        PickBuffS(Pick) = 1
    ElseIf Flag = 0 Then
        Erase PickBuffS
        Pick = Int(Rnd * Size) + 1
        PickBuffS(Pick) = 1
    End If
End Sub


Private Sub subAddProbeRecord(lTrial As Long)
    
    datProbe.Recordset.AddNew
    datProbe.Recordset("TrialNum") = lTrial
    datProbe.Recordset("ViewTime") = datSessionInfo.Recordset("DefaultViewTime")
    datProbe.Recordset("Delay") = datSessionInfo.Recordset("DefaultProbeDelay")
    datProbe.Recordset.Update

End Sub

Private Sub subAddStimuliRecords(lTrial As Long, ByVal randDRO As Boolean)
    
    Dim iPosition As Integer
        
    For iPosition = 1 To datSessionInfo.Recordset("DefaultMaxListLength")
        
        DoEvents
        If bLoadingRecords = False Then
            Exit Sub
        End If
        
        datStimuli.Recordset.AddNew
        datStimuli.Recordset("TrialNum") = lTrial
        datStimuli.Recordset("Position") = iPosition
        datStimuli.Recordset("ISI") = datSessionInfo.Recordset("DefaultISI")
        
        If randDRO = True Then
            'choose randomly between 1-5s view times
            Dim low As Integer
            Dim high As Integer
            low = 1000
            high = 5000
            
            Dim choice As Integer
            'generate 0-8
            choice = Int(Rnd() * 9)
            
            Debug.Print choice
            
            datStimuli.Recordset("ViewTime") = low + 500 * choice
        Else
            datStimuli.Recordset("ViewTime") = 1000
        End If
        
        datStimuli.Recordset.Update
    Next iPosition
    
End Sub

Private Sub subAddTrialRecord(lTrial As Long, ByVal useDRO As Boolean)
    datTrialInfo.Recordset.AddNew
    datTrialInfo.Recordset("TrialNum") = lTrial
    'If useDRO Then
        'datTrialInfo.Recordset("TrialType") = "DRO"
    'Else
         datTrialInfo.Recordset("TrialType") = "NO_DRO"
    'End If
    datTrialInfo.Recordset("ListLength") = datSessionInfo.Recordset("DefaultMaxListLength")
    datTrialInfo.Recordset("CorrectPosition") = 0
    
    
    'ITI is U[min, max] unless min = max, then just assign as min
    Dim iti, min, max As Integer
    min = datSessionInfo.Recordset("ProbeTime")
    max = datSessionInfo.Recordset("ProbeAbortLTO")
    If min = max Then
        iti = min
    Else
        Dim range, choice As Integer
        range = (max - min) + 1
        choice = Int(Rnd() * range)
        
        iti = min + choice
    End If
    datTrialInfo.Recordset("ITI") = iti
    
    datTrialInfo.Recordset("Reinforce") = True
    datTrialInfo.Recordset("Apparatus") = "Pellets"
    If datSessionInfo.Recordset("NoTouchProbes") = True Then datTrialInfo.Recordset("NoTouchTrial") = True
    datTrialInfo.Recordset.Update
End Sub


Private Sub subAutoAssignPosition()
    Dim i, display_size, num_positions As Integer
    display_size = 6
    num_positions = 6
    Dim iTrials As Integer
    iTrials = datSessionInfo.Recordset("NumOfTrials")
    datTrialInfo.Recordset.MoveFirst
    For i = 1 To iTrials
            datTrialInfo.Recordset.Edit
        Dim pos(6) As Integer
        pos(1) = Int((num_positions * Rnd) + 1)
        
            Dim c, cc As Integer
            Dim go_on As Boolean
            
            For c = 2 To display_size
                Do
                    pos(c) = Int((num_positions * Rnd) + 1)
                    cc = 1
                    Do
                        go_on = pos(cc) <> pos(c)
                        cc = cc + 1
                    Loop While cc < c And go_on = True
                Loop While Not go_on
            Next c
        
        
        'add all the positions to the trialInfo table
        For c = 1 To display_size
            datTrialInfo.Recordset("POS" & c) = pos(c)
        Next c
        
        'because all locations are randomly assigned, it doesn't matter which
        'position is the one that changed, let's go with pos(1)
        'in fact, subDisplaySecond is set to bail on POSChange <> Pos1
        datTrialInfo.Recordset("POSchange") = pos(1)
        datTrialInfo.Recordset.Update
        datTrialInfo.Recordset.MoveNext
     Next i
End Sub


Private Sub subAssignApparatus()

    Dim iTrials As Integer
    Dim iSameTrials As Integer
'    Dim iDiffTrials As Integer
    Dim iRecord As Integer
    Dim i As Integer
    Dim iSameNoReinforce As Integer
'    Dim iDiffNoReinforce As Integer
    Dim iSameJuice As Integer
    Dim iSamePellets As Integer
'    Dim iDiffJuice As Integer
'    Dim iDiffPellets As Integer
    Dim sngJuiceProb As Single
    Dim sngPelletProb As Single

    
    iTrials = datSessionInfo.Recordset("NumOfTrials")
    
    'just use pellet probability to control everything
    iSameTrials = iTrials ' * 0.5  'datSessionInfo.Recordset("SameProbability") \ 100
    iSamePellets = iSameTrials * (datSessionInfo.Recordset("PelletProbability") / 100)
    iSameJuice = iSameTrials - iSamePellets
    
'    iDiffTrials = iTrials - iSameTrials
'    iDiffPellets = iDiffTrials * (datSessionInfo.Recordset("PelletProbability") / 100)
'    iDiffJuice = iDiffTrials - iDiffPellets
    
    Randomize
    frameStatus.ZOrder
    frameStatus.Visible = True
    
    For i = 1 To iSameJuice
    
        frameStatus.Caption = "Assigning Same Juice:  " & i & " of " & iSameJuice
        Refresh
        
        Do
            DoEvents
            If bLoadingRecords = False Then
                Exit Sub
            End If

            iRecord = Int((iTrials * Rnd) + 1)
            datTrialInfo.Recordset.FindFirst "TrialNum = " & iRecord
        Loop While datTrialInfo.Recordset("Apparatus") = "Juice" Or datTrialInfo.Recordset("Apparatus") = "Pellets" Or datTrialInfo.Recordset("TrialType") = "Diff" Or datTrialInfo.Recordset("Reinforce") = False
        
        datTrialInfo.Recordset.Edit
        datTrialInfo.Recordset("Apparatus") = "Juice"
        datTrialInfo.Recordset.Update
        
    Next i
    
    
    For i = 1 To iSamePellets
    
        frameStatus.Caption = "Assigning Same Juice:  " & i & " of " & iSamePellets
        Refresh
        
        Do
            DoEvents
            If bLoadingRecords = False Then
                Exit Sub
            End If
            iRecord = Int((iTrials * Rnd) + 1)
            datTrialInfo.Recordset.FindFirst "TrialNum = " & iRecord
        Loop While datTrialInfo.Recordset("Apparatus") = "Juice" Or datTrialInfo.Recordset("Apparatus") = "Pellets" Or datTrialInfo.Recordset("TrialType") = "Diff" Or datTrialInfo.Recordset("Reinforce") = False
        
        datTrialInfo.Recordset.Edit
        datTrialInfo.Recordset("Apparatus") = "Pellets"
        datTrialInfo.Recordset.Update
        
    Next i
    
    
'    For i = 1 To iDiffJuice
    
'        frameStatus.Caption = "Assigning Diff Juice:  " & i & " of " & iDiffJuice
'        Refresh
'
'        Do
'            DoEvents
'            If bLoadingRecords = False Then
'                Exit Sub
'            End If
'            iRecord = Int((iTrials * Rnd) + 1)
'            datTrialInfo.Recordset.FindFirst "TrialNum = " & iRecord
'        Loop While datTrialInfo.Recordset("Apparatus") = "Juice" Or datTrialInfo.Recordset("TrialType") = "Same" Or datTrialInfo.Recordset("Reinforce") = False
'
'        datTrialInfo.Recordset.Edit
'        datTrialInfo.Recordset("Apparatus") = "Juice"
'        datTrialInfo.Recordset.Update
'
'    Next i
'
'
'
'    For i = 1 To iDiffPellets
'
'        frameStatus.Caption = "Assigning Diff Pellets:  " & i & " of " & iDiffPellets
'        Refresh
'
'        Do
'            DoEvents
'            If bLoadingRecords = False Then
'                Exit Sub
'            End If
'            iRecord = Int((iTrials * Rnd) + 1)
'            datTrialInfo.Recordset.FindFirst "TrialNum = " & iRecord
'        Loop While datTrialInfo.Recordset("Apparatus") = "Juice" Or datTrialInfo.Recordset("Apparatus") = "Pellets" Or datTrialInfo.Recordset("TrialType") = "Same" Or datTrialInfo.Recordset("Reinforce") = False
'
'        datTrialInfo.Recordset.Edit
'        datTrialInfo.Recordset("Apparatus") = "Pellets"
'        datTrialInfo.Recordset.Update
'
'    Next i
    
    frameStatus.Visible = False
    
End Sub

Private Sub subAssignReinforce()

    Dim iTrials As Integer
    Dim iSameTrials As Integer
    Dim iDiffTrials As Integer
    Dim iReinfTrials As Integer
    Dim iRecord As Integer
    Dim i As Integer
    Dim iSameReinforce As Integer
    Dim iDiffReinforce As Integer
    
    iTrials = datSessionInfo.Recordset("NumOfTrials")
    iSameTrials = iTrials * datSessionInfo.Recordset("SameProbability") \ 100
    iDiffTrials = iTrials - iSameTrials
    
    iReinfTrials = iTrials '* datSessionInfo.Recordset("ReinfProbability") / 100
    iSameReinforce = iSameTrials '* datSessionInfo.Recordset("ReinfProbability") \ 100
    iDiffReinforce = iReinfTrials - iSameReinforce
    
    Randomize
    frameStatus.ZOrder
    frameStatus.Visible = True
    
    For i = 1 To iSameReinforce
    
        frameStatus.Caption = "Assigning Same Reinforce:  " & i & " of " & iSameReinforce
        Refresh
        
        Do
            DoEvents
            If bLoadingRecords = False Then
                Exit Sub
            End If
            iRecord = Int((iTrials * Rnd) + 1)
            datTrialInfo.Recordset.FindFirst "TrialNum = " & iRecord
        Loop While datTrialInfo.Recordset("TrialType") = "Diff" Or datTrialInfo.Recordset("Reinforce") = True
        
        datTrialInfo.Recordset.Edit
        datTrialInfo.Recordset("Reinforce") = True
        datTrialInfo.Recordset("Apparatus") = "None"
        datTrialInfo.Recordset.Update
        
    Next i
    
    For i = 1 To iDiffReinforce
    
        frameStatus.Caption = "Assigning Diff Reinforce:  " & i & " of " & iDiffReinforce
        Refresh
        
        Do
            DoEvents
            If bLoadingRecords = False Then
                Exit Sub
            End If
            iRecord = Int((iTrials * Rnd) + 1)
            datTrialInfo.Recordset.FindFirst "TrialNum = " & iRecord
        Loop While datTrialInfo.Recordset("TrialType") = "Same" Or datTrialInfo.Recordset("Reinforce") = True
        
        datTrialInfo.Recordset.Edit
        datTrialInfo.Recordset("Reinforce") = True
        datTrialInfo.Recordset("Apparatus") = "None"
        datTrialInfo.Recordset.Update
        
    Next i
    
    frameStatus.Visible = False
    
End Sub
Private Sub subAssignConfiguration()
    Dim iTrials As Integer
    Dim i, t, tally As Integer
    Dim Pick As Integer
    Dim Flag As Boolean
    Dim range As Integer
    Dim Tset, Sset, Dset As Integer
    
    
    Randomize
    iTrials = datSessionInfo.Recordset("NumOfTrials")
    datTrialInfo.Recordset.MoveFirst
    frameStatus.ZOrder
    frameStatus.Visible = True
    
    Tset = 60
    Sset = 30
    Dset = 30
    
    Open "Stimuli.txt" For Input As #1
        Do While Not EOF(1) ' Loop until end of file.
            Input #1, i, t
            If i < (Sset + 1) Then
                PickBuffS(i) = t
                Debug.Print i, PickBuffS(i)  ' Print data to Debug window.
            Else
                PickBuffD(i - Sset) = t
                Debug.Print i, PickBuffD(i - Sset)
            End If
        Loop
            Debug.Print tally
    Close #1    ' Close file.
    
    For i = 1 To iTrials
        frameStatus.Caption = "Assigning Configuration:  " & i & " of " & iTrials
        Refresh
        datTrialInfo.Recordset.Edit
        Select Case datTrialInfo.Recordset("TrialType")
                Case "Left"
                        PickWoReplaceS 30, Pick
                        datTrialInfo.Recordset("Configuration") = Pick
                    
                Case "Right"
                        PickWoReplaceD Dset, Pick
                        datTrialInfo.Recordset("Configuration") = Pick + Sset
                Case Else
                    Stop
        End Select
        
        datTrialInfo.Recordset("Reinforce") = True
        datTrialInfo.Recordset("Apparatus") = "Pellets"
        datTrialInfo.Recordset.Update
        datTrialInfo.Recordset.MoveNext
     Next i
     
     Open "Stimuli.txt" For Output Shared As #1
            For i = 1 To Tset
                If i < (Sset + 1) Then
                    Write #1, i; PickBuffS(i)
                Else
                    Write #1, i; PickBuffD(i - Sset)
                End If
            Next i
        Close #1
    frameStatus.Visible = False

End Sub
Private Sub subAssignPosition1(ps As Integer)

Select Case ps
            Case 1
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 2
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 2
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 2
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 3
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 3
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 4
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 3
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 5
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 4
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 6
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 4
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 7
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 5
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 8
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 5
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 9
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 6
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 10
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 6
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 11
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 7
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 12
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 7
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 13
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 8
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 14
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 8
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 15
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 9
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 16
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 9
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 17
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 10
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 18
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 10
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 19
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 20
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 21
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 22
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 23
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 24
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 25
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 26
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 27
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 28
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 29
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 30
                datTrialInfo.Recordset("POS1") = 1
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 31
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 3
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 32
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 3
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 33
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 4
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 34
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 4
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 35
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 5
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 36
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 5
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 37
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 6
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 38
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 6
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 39
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 7
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 40
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 7
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 41
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 8
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 42
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 8
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 43
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 9
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 44
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 9
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 45
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 10
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 46
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 10
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 47
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 48
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 49
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 50
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 51
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 52
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 53
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 54
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 55
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 56
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 57
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 58
                datTrialInfo.Recordset("POS1") = 2
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 59
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 4
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 60
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 4
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 61
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 5
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 62
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 5
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 63
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 6
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 64
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 6
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 65
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 7
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 66
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 7
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 67
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 8
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 68
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 8
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 69
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 9
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 70
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 9
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 71
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 10
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 72
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 10
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 73
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 74
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 75
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 76
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 77
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 78
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 79
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 80
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 81
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 82
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 83
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 84
                datTrialInfo.Recordset("POS1") = 3
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 85
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 5
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 86
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 5
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 87
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 6
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 88
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 6
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 89
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 7
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 90
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 7
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 91
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 8
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 92
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 8
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 93
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 9
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 94
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 9
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 95
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 10
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 96
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 10
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 97
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 98
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 99
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 100
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 101
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 102
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 103
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 104
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 105
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 106
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 107
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 108
                datTrialInfo.Recordset("POS1") = 4
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 109
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 6
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 110
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 6
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 111
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 7
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 112
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 7
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 113
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 8
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 114
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 8
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 115
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 9
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 116
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 9
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 117
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 10
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 118
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 10
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 119
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 120
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case Else
                Stop
        End Select


        datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")

End Sub
Private Sub subAssignPosition2(ps As Integer)
                
        Select Case ps
            Case 121
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 122
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 123
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 124
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 125
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 126
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 127
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 128
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 129
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 130
                datTrialInfo.Recordset("POS1") = 5
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 131
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 7
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 132
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 7
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 133
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 8
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 134
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 8
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 135
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 9
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 136
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 9
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 137
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 10
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 138
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 10
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 139
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 140
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 141
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 142
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 143
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 144
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 145
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 146
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 147
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 148
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 149
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 150
                datTrialInfo.Recordset("POS1") = 6
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 151
                datTrialInfo.Recordset("POS1") = 7
                datTrialInfo.Recordset("POS2") = 8
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 152
                datTrialInfo.Recordset("POS1") = 7
                datTrialInfo.Recordset("POS2") = 8
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 153
                datTrialInfo.Recordset("POS1") = 7
                datTrialInfo.Recordset("POS2") = 9
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 154
                datTrialInfo.Recordset("POS1") = 7
                datTrialInfo.Recordset("POS2") = 9
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 155
                datTrialInfo.Recordset("POS1") = 7
                datTrialInfo.Recordset("POS2") = 10
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 156
                datTrialInfo.Recordset("POS1") = 7
                datTrialInfo.Recordset("POS2") = 10
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 157
                datTrialInfo.Recordset("POS1") = 7
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 158
                datTrialInfo.Recordset("POS1") = 7
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 159
                datTrialInfo.Recordset("POS1") = 7
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 160
                datTrialInfo.Recordset("POS1") = 7
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 161
                datTrialInfo.Recordset("POS1") = 7
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 162
                datTrialInfo.Recordset("POS1") = 7
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 163
                datTrialInfo.Recordset("POS1") = 7
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 164
                datTrialInfo.Recordset("POS1") = 7
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 165
                datTrialInfo.Recordset("POS1") = 7
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 166
                datTrialInfo.Recordset("POS1") = 7
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 167
                datTrialInfo.Recordset("POS1") = 7
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 168
                datTrialInfo.Recordset("POS1") = 7
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 169
                datTrialInfo.Recordset("POS1") = 8
                datTrialInfo.Recordset("POS2") = 9
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 170
                datTrialInfo.Recordset("POS1") = 8
                datTrialInfo.Recordset("POS2") = 9
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 171
                datTrialInfo.Recordset("POS1") = 8
                datTrialInfo.Recordset("POS2") = 10
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 172
                datTrialInfo.Recordset("POS1") = 8
                datTrialInfo.Recordset("POS2") = 10
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 173
                datTrialInfo.Recordset("POS1") = 8
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 174
                datTrialInfo.Recordset("POS1") = 8
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 175
                datTrialInfo.Recordset("POS1") = 8
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 176
                datTrialInfo.Recordset("POS1") = 8
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 177
                datTrialInfo.Recordset("POS1") = 8
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 178
                datTrialInfo.Recordset("POS1") = 8
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 179
                datTrialInfo.Recordset("POS1") = 8
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 180
                datTrialInfo.Recordset("POS1") = 8
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 181
                datTrialInfo.Recordset("POS1") = 8
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 182
                datTrialInfo.Recordset("POS1") = 8
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 183
                datTrialInfo.Recordset("POS1") = 8
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 184
                datTrialInfo.Recordset("POS1") = 8
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 185
                datTrialInfo.Recordset("POS1") = 9
                datTrialInfo.Recordset("POS2") = 10
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 186
                datTrialInfo.Recordset("POS1") = 9
                datTrialInfo.Recordset("POS2") = 10
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 187
                datTrialInfo.Recordset("POS1") = 9
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 188
                datTrialInfo.Recordset("POS1") = 9
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 189
                datTrialInfo.Recordset("POS1") = 9
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 190
                datTrialInfo.Recordset("POS1") = 9
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 191
                datTrialInfo.Recordset("POS1") = 9
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 192
                datTrialInfo.Recordset("POS1") = 9
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 193
                datTrialInfo.Recordset("POS1") = 9
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 194
                datTrialInfo.Recordset("POS1") = 9
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 195
                datTrialInfo.Recordset("POS1") = 9
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 196
                datTrialInfo.Recordset("POS1") = 9
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 197
                datTrialInfo.Recordset("POS1") = 9
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 198
                datTrialInfo.Recordset("POS1") = 9
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 199
                datTrialInfo.Recordset("POS1") = 10
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 200
                datTrialInfo.Recordset("POS1") = 10
                datTrialInfo.Recordset("POS2") = 11
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 201
                datTrialInfo.Recordset("POS1") = 10
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 202
                datTrialInfo.Recordset("POS1") = 10
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 203
                datTrialInfo.Recordset("POS1") = 10
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 204
                datTrialInfo.Recordset("POS1") = 10
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 205
                datTrialInfo.Recordset("POS1") = 10
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 206
                datTrialInfo.Recordset("POS1") = 10
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 207
                datTrialInfo.Recordset("POS1") = 10
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 208
                datTrialInfo.Recordset("POS1") = 10
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 209
                datTrialInfo.Recordset("POS1") = 10
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 210
                datTrialInfo.Recordset("POS1") = 10
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 211
                datTrialInfo.Recordset("POS1") = 11
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 212
                datTrialInfo.Recordset("POS1") = 11
                datTrialInfo.Recordset("POS2") = 12
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 213
                datTrialInfo.Recordset("POS1") = 11
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 214
                datTrialInfo.Recordset("POS1") = 11
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 215
                datTrialInfo.Recordset("POS1") = 11
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 216
                datTrialInfo.Recordset("POS1") = 11
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 217
                datTrialInfo.Recordset("POS1") = 11
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 218
                datTrialInfo.Recordset("POS1") = 11
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 219
                datTrialInfo.Recordset("POS1") = 11
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 220
                datTrialInfo.Recordset("POS1") = 11
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 221
                datTrialInfo.Recordset("POS1") = 12
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 222
                datTrialInfo.Recordset("POS1") = 12
                datTrialInfo.Recordset("POS2") = 13
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 223
                datTrialInfo.Recordset("POS1") = 12
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 224
                datTrialInfo.Recordset("POS1") = 12
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 225
                datTrialInfo.Recordset("POS1") = 12
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 226
                datTrialInfo.Recordset("POS1") = 12
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 227
                datTrialInfo.Recordset("POS1") = 12
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 228
                datTrialInfo.Recordset("POS1") = 12
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 229
                datTrialInfo.Recordset("POS1") = 13
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 230
                datTrialInfo.Recordset("POS1") = 13
                datTrialInfo.Recordset("POS2") = 14
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 231
                datTrialInfo.Recordset("POS1") = 13
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 232
                datTrialInfo.Recordset("POS1") = 13
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 233
                datTrialInfo.Recordset("POS1") = 13
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 234
                datTrialInfo.Recordset("POS1") = 13
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 235
                datTrialInfo.Recordset("POS1") = 14
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 236
                datTrialInfo.Recordset("POS1") = 14
                datTrialInfo.Recordset("POS2") = 15
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 237
                datTrialInfo.Recordset("POS1") = 14
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 238
                datTrialInfo.Recordset("POS1") = 14
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case 239
                datTrialInfo.Recordset("POS1") = 15
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
            Case 240
                datTrialInfo.Recordset("POS1") = 15
                datTrialInfo.Recordset("POS2") = 16
                datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS2")
            Case Else
                Stop
        End Select

datTrialInfo.Recordset("POSchange") = datTrialInfo.Recordset("POS1")
End Sub

Private Sub subAssignConfigurationOld()
    Dim iTrials As Integer
    Dim i, q, tally, tally1, Sets As Integer
    Dim Pick As Integer
        
    Randomize
    iTrials = datSessionInfo.Recordset("NumOfTrials")
    datTrialInfo.Recordset.MoveFirst
    frameStatus.ZOrder
    frameStatus.Visible = True
           
    Open "Location.txt" For Input As #1
        Do While Not EOF(1) ' Loop until end of file.
            Input #1, i, q
            PickBuffS(i) = q
            Debug.Print i, PickBuffS(i)  ' Print data to Debug window.
            If PickBuffS(i) = 1 Then tally = tally + 1
        Loop
            Debug.Print tally
    Close #1    ' Close file.
    If tally >= 240 Then Erase PickBuffS
   
    Open "Stimuli.txt" For Input As #1
        Do While Not EOF(1) ' Loop until end of file.
            Input #1, i, q
            PickBuffD(i) = q
            Debug.Print i, PickBuffD(i)  ' Print data to Debug window.
            If PickBuffD(i) = 1 Then tally1 = tally1 + 1
        Loop
            Debug.Print tally1
    Close #1    ' Close file.
    If tally1 >= 168 Then Erase PickBuffD

        Sets = datSessionInfo.Recordset("SetSize")
        If Sets = 0 Then
            ErrorBox.Text = " Set Size = 0. Restart"
            ErrorBox.Visible = True
            Exit Sub
        End If
        If Sets = 4 Then Erase PickBuffD
        
    For i = 1 To iTrials
        frameStatus.Caption = "Assigning Configuration:  " & i & " of " & iTrials
        Refresh
        datTrialInfo.Recordset.Edit
        
        If Sets = 4 Then
            Pick = Int((12 * Rnd) + 1) ' 12 stimulus types
            PickWoReplaceD 12, Pick
            datTrialInfo.Recordset("Configuration") = Pick
        ElseIf Sets = 8 Then
            Pick = Int((168 * Rnd) + 1) ' 168 stimulus types
            PickWoReplaceD 168, Pick
            datTrialInfo.Recordset("Configuration") = Pick
        Else
            Stop
        End If
                        
        Pick = Int((240 * Rnd) + 1)
        PickWoReplaceS 240, Pick
                                             
        Select Case Pick
            Case 1 To 120
                subAssignPosition1 (Pick)
            Case 121 To 240
                subAssignPosition2 (Pick)
        End Select
        
        datTrialInfo.Recordset("Reinforce") = True
        datTrialInfo.Recordset("Apparatus") = "Pellets"
        datTrialInfo.Recordset.Update
        datTrialInfo.Recordset.MoveNext
     Next i
     tally = 0
     Open "Location.txt" For Output Shared As #1
         For i = 1 To 240
            Print #1, i; PickBuffS(i)
            If PickBuffS(i) = 1 Then tally = tally + 1
        Next i
        Debug.Print tally
     Close #1
     tally1 = 0
     Open "Stimuli.txt" For Output Shared As #1
        For i = 1 To 168
            Print #1, i; PickBuffD(i)
            If PickBuffD(i) = 1 Then tally1 = tally1 + 1
        Next i
        Debug.Print tally1
     Close #1
    
     frameStatus.Visible = False

End Sub
Private Sub subAssignTypesLeftRight()

    Dim iTrials As Integer
    Dim iSameTrials As Integer  'SAME = LEFT
    Dim iDiffTrials As Integer  'DIFFERENT = RIGHT
    Dim iRecord As Integer
    Dim i As Integer

    iTrials = datSessionInfo.Recordset("NumOfTrials")
    iSameTrials = iTrials * datSessionInfo.Recordset("SameProbability") \ 100
    'adjust for adding trials to an existing session
    iDiffTrials = iTrials - iSameTrials - CInt(lblDiffsCount.Caption)
    
    Randomize
    frameStatus.ZOrder
    frameStatus.Visible = True
    For i = 1 To iDiffTrials
    
        frameStatus.Caption = "Assigning Rights:  " & i & " of " & iDiffTrials
        Refresh
        
        Do
            DoEvents
            If bLoadingRecords = False Then
                Exit Sub
            End If
            iRecord = Int((iTrials * Rnd) + 1)
            datTrialInfo.Recordset.FindFirst "TrialNum = " & iRecord
        Loop While datTrialInfo.Recordset("TrialType") = "Right"
        
        datTrialInfo.Recordset.Edit
        datTrialInfo.Recordset("TrialType") = "Right"
        datTrialInfo.Recordset.Update
        
    Next i
    
    frameStatus.Visible = False
    

End Sub

Private Sub subAssignTypesSameDifferent()

    Dim iTrials As Integer
    Dim iSameTrials As Integer
    Dim iDiffTrials As Integer
    Dim iRecord As Integer
    Dim i As Integer

    iTrials = datSessionInfo.Recordset("NumOfTrials")
    iSameTrials = iTrials * datSessionInfo.Recordset("SameProbability") \ 100
    'adjust for adding trials to an existing session
    iDiffTrials = iTrials - iSameTrials - CInt(lblDiffsCount.Caption)
    
    Randomize
    frameStatus.ZOrder
    frameStatus.Visible = True
    For i = 1 To iDiffTrials
    
        frameStatus.Caption = "Assigning Diffs:  " & i & " of " & iDiffTrials
        Refresh
        
        Do
            DoEvents
            If bLoadingRecords = False Then
                Exit Sub
            End If
            iRecord = Int((iTrials * Rnd) + 1)
            datTrialInfo.Recordset.FindFirst "TrialNum = " & iRecord
        Loop While datTrialInfo.Recordset("TrialType") = "Diff"
        
        datTrialInfo.Recordset.Edit
        datTrialInfo.Recordset("TrialType") = "Diff"
        datTrialInfo.Recordset.Update
        
    Next i
    
    frameStatus.Visible = False
    

End Sub

Private Sub subCalcApparatusPercentages()

    If (CInt(lblJuiceCount.Caption) + CInt(lblPelletsCount.Caption)) <> 0 Then
        lblJuiceActual.Caption = CInt(lblJuiceCount.Caption) * 100 \ (CInt(lblJuiceCount.Caption) + CInt(lblPelletsCount.Caption))
        lblPelletsActual.Caption = CInt(lblPelletsCount.Caption) * 100 \ (CInt(lblJuiceCount.Caption) + CInt(lblPelletsCount.Caption))
    Else
        lblJuiceActual.Caption = 0
        lblPelletsActual.Caption = 0
    End If
    
End Sub

Private Sub subCalcReinforcePercentages()

    If CInt(lblNumTrials.Caption) <> 0 Then
        lblReinforceActual.Caption = CInt(lblReinforceCount.Caption) * 100 \ CInt(lblNumTrials.Caption)
    Else
        lblReinforceActual.Caption = 0
    End If
    
End Sub


Private Sub subCalcTypePercentages()

    If CInt(lblNumTrials.Caption) <> 0 Then
        lblDiffsActual.Caption = CInt(lblDiffsCount.Caption) * 100 \ CInt(lblNumTrials.Caption)
        lblSamesActual.Caption = CInt(lblSamesCount.Caption) * 100 \ CInt(lblNumTrials.Caption)
    Else
        lblDiffsActual.Caption = 0
        lblSamesActual.Caption = 0
    End If
    
End Sub

Private Sub subClearApparatus()
    
    datTrialInfo.Recordset.MoveFirst
    
    Do While Not datTrialInfo.Recordset.EOF
    
        DoEvents
        
        datTrialInfo.Recordset.Edit
        datTrialInfo.Recordset("Apparatus") = "None"
        datTrialInfo.Recordset.Update
        
        datTrialInfo.Recordset.MoveNext
        
    Loop
    
End Sub

Private Sub subCountDiffs()

    datRecordCounter.DatabaseName = gsProfile
    datRecordCounter.RecordSource = "SELECT COUNT(*) as Quantity FROM TrialInfo WHERE TrialType = 'Right'"
    datRecordCounter.Refresh

    lblDiffsCount.Caption = datRecordCounter.Recordset("Quantity")
    
End Sub

Private Sub subCountJuice()
    
    datRecordCounter.DatabaseName = gsProfile
    datRecordCounter.RecordSource = "SELECT COUNT(*) as Quantity FROM TrialInfo WHERE Apparatus = 'Juice'"
    datRecordCounter.Refresh

    lblJuiceCount.Caption = datRecordCounter.Recordset("Quantity")
    
End Sub

Private Sub subCountPellets()
    
    datRecordCounter.DatabaseName = gsProfile
    datRecordCounter.RecordSource = "SELECT COUNT(*) as Quantity FROM TrialInfo WHERE Apparatus = 'Pellets'"
    datRecordCounter.Refresh

    lblPelletsCount.Caption = datRecordCounter.Recordset("Quantity")

End Sub

Private Sub subCountReinforces()
    
    datRecordCounter.DatabaseName = gsProfile
    datRecordCounter.RecordSource = "SELECT COUNT(*) as Quantity FROM TrialInfo WHERE Reinforce = TRUE"
    datRecordCounter.Refresh

    lblReinforceCount.Caption = datRecordCounter.Recordset("Quantity")
    
End Sub

Private Sub subCountSames()
    
    datRecordCounter.DatabaseName = gsProfile
    datRecordCounter.RecordSource = "SELECT COUNT(*) as Quantity FROM TrialInfo WHERE TrialType = 'Left'"
    datRecordCounter.Refresh

    lblSamesCount.Caption = datRecordCounter.Recordset("Quantity")
    
End Sub

Private Sub subLoadTrialRecords()

    Dim lNumRecords As Long
    Dim lNumTrials As Long
    Dim lTrial As Long
    Dim i As Long
    
    
    'turn off status
    gridStimuli.Visible = False
    gridProbe.Visible = False
    gridTrialInfo.Visible = False
    frameStatus.ZOrder
    frameStatus.Visible = True
    
    'check for correct number of records
    If datTrialInfo.Recordset.RecordCount <> 0 Then
        datTrialInfo.Recordset.MoveLast
        lTrial = datTrialInfo.Recordset("TrialNum")
    Else
        lTrial = 0
    End If
    
    lNumRecords = datTrialInfo.Recordset.RecordCount
    lNumTrials = datSessionInfo.Recordset("NumOfTrials")
    
    Dim useDRO(1 To 192) As Integer
    Dim d As Integer

    For d = 1 To 192
        useDRO(d) = False
    Next d

    For d = 1 To 96
        Do
            Dim choice As Integer
            choice = Int(Rnd * 192 + 1)
        Loop While (useDRO(choice) = True)

        useDRO(choice) = True
    Next d
    
    If lNumRecords < lNumTrials Then
    
        For i = 1 To (lNumTrials - lNumRecords)
            DoEvents
            If bLoadingRecords = False Then
                Exit Sub
            End If
            
            lTrial = lTrial + 1
            frameStatus.Caption = "Loading Records:  " & i & " of " & (lNumTrials - lNumRecords)
            subAddTrialRecord lTrial, useDRO(lTrial)
            
            If datSessionInfo.Recordset("DefaultVariableListLength") = False Then
                subAddStimuliRecords lTrial, useDRO(lTrial)
            End If
            subAddProbeRecord lTrial
        Next i
        datTrialInfo.Refresh
        datTrialInfo.Recordset.MoveFirst
        
        datStimuli.RecordSource = "SELECT * FROM StimuliInfo WHERE TrialNum = 1"
        datStimuli.Refresh
        datProbe.RecordSource = "SELECT * FROM ProbeInfo WHERE TrialNum = 1"
        datProbe.Refresh
    
        'note that we are assigning variable probe delays
        subAssignVariableProbeDelays
        
        'subAssignConfigurationOld
        subAutoAssignPosition
        subAssignDisplaySize 'uses the configuration column
        subAssignStim
   
        'choose apparatus
        subClearApparatus
        subAssignApparatus
    
    ElseIf lNumRecords > lNumTrials Then
        
        For i = 1 To (lNumRecords - lNumTrials)
            DoEvents
            If bLoadingRecords = False Then
                Exit Sub
            End If

            'delete Trial record
            datTrialInfo.Recordset.Delete
            datTrialInfo.Recordset.MoveLast
        Next i
    
    End If
    
    datTrialInfo.Refresh
    
    'turn off status
    frameStatus.Visible = False
    gridStimuli.Visible = True
    gridProbe.Visible = True
    gridTrialInfo.Visible = True


End Sub


Private Sub subAssignDisplaySize()
    Dim lNumTrials As Integer
    lNumTrials = datSessionInfo.Recordset("NumOfTrials")
    
    'we need an array that we can store with various display sizes
    'zero out the trial types array
    Dim trialTypes(192) As Integer
    
    Dim tt As Integer
    For tt = 1 To lNumTrials
        trialTypes(tt) = -1
    Next tt

    'names are created as array to facilitate storage later
    'for ds2-12 we need 11 spaces
    Dim trialTypeNames(11) As Integer
    
    'based on the values stored in ds(i) we need to set i
    'display sizes on a trial-by-trial basis
    Dim displaySizes(11) As Integer
    For tt = 1 To 11
        displaySizes(tt) = Int(datView.Recordset("ds" & (tt + 1)))
        trialTypeNames(tt) = (tt + 1)
    Next tt
    
    Dim randDS As Boolean
    randDS = False
    If (datSessionInfo.Recordset("TrialUnique") = 0) Then
        randDS = True
    End If
    
    'now that we have the number of each display size type,
    'we just loop through each size and assign a random trial to that type,
    'making sure that the type of that trial was "NONE"
    'the trialTypes array will be used when we actually load trials
    Dim i As Integer
    i = 1
    For tt = 1 To 11
        Dim ds As Integer
        For ds = 1 To displaySizes(tt)
            Debug.Print "building display size " & ds & ": " & displaySizes(tt)
            
            If (randDS) Then

                'pick a random trial to set to the current display size
                Do
                    Dim randTrial As Integer
                    randTrial = Int((Rnd * lNumTrials) + 1)
                
                'we keep looping if the trial type wasn't "NONE"
                'this could be replaced by a better sample w/o replacement alg.
                Loop While trialTypes(randTrial) <> -1
                
                'because the trialType at randTrial = "NONE" we can safely
                'set it to the appropriate display size, using the names array
                'Note that we are using the outermost loop index variable `tt`
                trialTypes(randTrial) = trialTypeNames(tt)
            Else
                trialTypes(i) = trialTypeNames(tt)
            End If
            i = i + 1
        Next ds
    Next tt
    Debug.Print "built trialTypes array"
    
    Dim c As Integer
    For c = 1 To lNumTrials
        datTrialInfo.Recordset.FindFirst "TrialNum = " & c
        datTrialInfo.Recordset.Edit
        datTrialInfo.Recordset("Configuration") = trialTypes(c)
        datTrialInfo.Recordset.Update

    Next c
End Sub

Private Sub subPrintSetup()
    
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
    
    bPrinting = True
    
    sField1 = Space(6)
    sField2 = sField1
    sField3 = sField1
    sField4 = sField1
    sField5 = sField1
    sField6 = sField1
    sField7 = sField1
    
        
    If datTrialInfo.Recordset.RecordCount <> 0 Then

        Printer.Font.name = "Courier"
        
        Printer.Print
        Printer.Print
        Printer.Print
        Printer.Print Tab(10); "SPR Setup Parameter: " & gsProfile
        
        Printer.Print
        Printer.Print
        Printer.Print
        
        RSet sField1 = lblReinforceCount
        RSet sField2 = lblReinforceActual & " %"
        RSet sField3 = lblReinforceDefault & " %"
        
        Printer.Print Tab(20); "  Count"; Tab(30); "Actual"; Tab(40); "Default"
        Printer.Print Tab(10); "Reinforce"; Tab(20); sField1; Tab(30); sField2; Tab(40); sField3
        
        Printer.Print
        Printer.Print
        Printer.Print
        
        RSet sField1 = lblJuiceCount
        RSet sField2 = lblJuiceActual & " %"
        RSet sField3 = lblJuiceDefault & " %"
        
        Printer.Print Tab(20); "  Count"; Tab(30); "Actual"; Tab(40); "Default"
        Printer.Print Tab(10); "Juice"; Tab(20); sField1; Tab(30); sField2; Tab(40); sField3
        
        RSet sField1 = lblPelletsCount
        RSet sField2 = lblPelletsActual & " %"
        RSet sField3 = lblPelletsDefault & " %"
        
        Printer.Print Tab(10); "Pellets"; Tab(20); sField1; Tab(30); sField2; Tab(40); sField3
        
        Printer.Print
        Printer.Print
        Printer.Print
        
        RSet sField1 = lblDiffsCount
        RSet sField2 = lblDiffsActual & " %"
        RSet sField3 = lblDiffsDefault & " %"
        
        Printer.Print Tab(20); "  Count"; Tab(30); "Actual"; Tab(40); "Default"
        Printer.Print Tab(10); "Diffs"; Tab(20); sField1; Tab(30); sField2; Tab(40); sField3
        
        RSet sField1 = lblSamesCount
        RSet sField2 = lblSamesActual & " %"
        RSet sField3 = lblSamesDefault & " %"
        
        Printer.Print Tab(10); "Sames"; Tab(20); sField1; Tab(30); sField2; Tab(40); sField3

        Printer.Print
        Printer.Print
        Printer.Print
        
        Printer.Print Tab(10); "Statistics"; Tab(30); "Position"; Tab(40); "  Count"
        Printer.Print
        
        datStats.Recordset.MoveFirst
        While Not datStats.Recordset.EOF
        
            DoEvents
            
            RSet sField1 = datStats.Recordset("CorrectPosition")
            RSet sField2 = datStats.Recordset("Quantity")
            Printer.Print Tab(30); sField1; _
                        Tab(40); sField2
        
            datStats.Recordset.MoveNext
            
        Wend
        
        Printer.Print
        Printer.Print
        Printer.Print
        
        Printer.Print Tab(10); "Total Trials = "; lblNumTrials

        datTrialInfo.Recordset.MoveFirst
        While Not datTrialInfo.Recordset.EOF
            
            DoEvents
            
            Printer.Print
            Printer.Print
            Printer.Print
            
            'heading
            Printer.Print Tab(10); "       "; Tab(20); "      "; Tab(30); "Correct"; Tab(40); "  List"
            Printer.Print Tab(10); "  Trial"; Tab(20); "  Type"; Tab(30); "Position"; Tab(40); "  Length"; Tab(50); "    ITI"; Tab(60); "Reinforce"; Tab(70); "Apparatus"
            Printer.Print
            
            RSet sField1 = lblTrialNum
            RSet sField2 = lblType
            RSet sField3 = intCorrectPosition
            RSet sField4 = intListLength
            RSet sField5 = intITI
            If chkReinforce <> 0 Then
                RSet sField6 = True
            Else
                RSet sField6 = False
            End If
            RSet sField7 = lblApparatus
            Printer.Print Tab(10); sField1; _
                        Tab(20); sField2; _
                        Tab(30); sField3; _
                        Tab(40); sField4; _
                        Tab(50); sField5; _
                        Tab(60); sField6; _
                        Tab(70); sField7
                        
                        
            Printer.Print
            Printer.Print
            Printer.Print
                        
            Printer.Print Tab(20); "Stimuli"; Tab(30); "Position"; Tab(40); "Filename"; Tab(60); "View Time"; Tab(70); "  ISI"
            Printer.Print
            
            datStimuli.Recordset.MoveFirst
            While Not datStimuli.Recordset.EOF
            
                DoEvents
                
                RSet sField1 = datStimuli.Recordset("Position")
                sField2 = datStimuli.Recordset("Filename")
                RSet sField3 = datStimuli.Recordset("ViewTime")
                RSet sField4 = datStimuli.Recordset("ISI")
                Printer.Print Tab(30); sField1; _
                            Tab(40); sField2; _
                            Tab(60); sField3; _
                            Tab(70); sField4
                            
                datStimuli.Recordset.MoveNext
                
            Wend
            
            Printer.Print
            Printer.Print
            Printer.Print
                        
            Printer.Print Tab(20); "Probe"; Tab(40); "Filename"; Tab(60); "View Time"; Tab(70); "  Delay"
            Printer.Print
            
            sField1 = datProbe.Recordset("Filename")
            RSet sField2 = datProbe.Recordset("ViewTime")
            RSet sField3 = datProbe.Recordset("Delay")
            Printer.Print Tab(40); sField1; _
                        Tab(60); sField2; _
                        Tab(70); sField3
            
            datTrialInfo.Recordset.MoveNext
            
        Wend
        
    End If
    
    Printer.EndDoc
    datTrialInfo.Recordset.MoveFirst
    bPrinting = False
    
End Sub

Private Sub chkReinforce_MouseUp(Button As Integer, Shift As Integer, X As Single, Y As Single)
    
    If bLoadingForm = False Then
        
        datTrialInfo.Recordset.Edit
        datTrialInfo.Recordset.Update
       
        subCountReinforces
        subCalcReinforcePercentages
        
        If chkReinforce = False Then
            datTrialInfo.Recordset.Edit
            datTrialInfo.Recordset("Apparatus") = "None"
            datTrialInfo.Recordset.Update
            subCountJuice
            subCountPellets
            subCalcApparatusPercentages
        End If
    
    End If

End Sub

Private Sub cmdDiff_Click()
    
    lblType.Caption = "Right"
    intCorrectPosition.Text = 0
    datTrialInfo.Recordset.Edit
    datTrialInfo.Recordset("TrialType") = "Right"
    datTrialInfo.Recordset.Update
    
    subCountSames
    subCountDiffs
    subCalculateStatistics
    subCalcTypePercentages
    
    Refresh
    
End Sub

Private Sub cmdJuice_Click()
    
    lblApparatus.Caption = "Juice"
    datTrialInfo.Recordset.Edit
    datTrialInfo.Recordset("Apparatus") = "Juice"
    datTrialInfo.Recordset.Update
   
    subCountJuice
    subCountPellets
    subCalcApparatusPercentages

End Sub

Private Sub cmdNextReinforceNoapparatus_Click()
    
    datTrialInfo.Recordset.FindNext "Reinforce = TRUE AND Apparatus = 'None'"

    If datTrialInfo.Recordset.NoMatch Then
        datTrialInfo.Recordset.FindFirst "Reinforce = TRUE AND Apparatus = 'None'"
    End If
    
End Sub

Private Sub cmdNextSame0Position_Click()

    datTrialInfo.Recordset.FindNext "TrialType = 'Same' AND CorrectPosition = 0"

    If datTrialInfo.Recordset.NoMatch Then
        datTrialInfo.Recordset.FindFirst "TrialType = 'Same' AND CorrectPosition = 0"
    End If
    
End Sub

Private Sub cmdPellets_Click()
    
    lblApparatus.Caption = "Pellets"
    datTrialInfo.Recordset.Edit
    datTrialInfo.Recordset("Apparatus") = "Pellets"
    datTrialInfo.Recordset.Update
   
    subCountJuice
    subCountPellets
    subCalcApparatusPercentages

End Sub


Private Sub cmdSame_Click()
    
    lblType.Caption = "Left"
    datTrialInfo.Recordset.Edit
    datTrialInfo.Recordset("TrialType") = "Left"
    datTrialInfo.Recordset.Update
   
    subCountSames
    subCountDiffs
    subCalculateStatistics
    subCalcTypePercentages
    
End Sub


Private Sub datTrialInfo_Reposition()
    
    If bLoadingForm = False Then
        
        'load the Stimuli and Probe tables into grids
        If datTrialInfo.Recordset.RecordCount > 0 And Not datTrialInfo.Recordset.EOF Then
            datStimuli.RecordSource = "SELECT * FROM StimuliInfo WHERE TrialNum = " & datTrialInfo.Recordset("TrialNum")
            datStimuli.Refresh
            datProbe.RecordSource = "SELECT * FROM ProbeInfo WHERE TrialNum = " & datTrialInfo.Recordset("TrialNum")
            datProbe.Refresh
        End If
        
        'make sure that the list length is not exceeded
        'intCorrectPosition.Max = intListLength.Text 'CLng(datTrialInfo.Recordset("ListLength"))
        
    End If
    
End Sub

Private Sub datTrialInfo_Validate(Action As Integer, Save As Integer)
    
    If bLoadingForm = False And bPrinting = False Then
        'make sure that new data is written to the table
        If datStimuli.Recordset.RecordCount <> 0 And Not datStimuli.Recordset.EOF Then
            datStimuli.Recordset.Edit
            datStimuli.Recordset.Update
        End If
        If datProbe.Recordset.RecordCount <> 0 And Not datProbe.Recordset.EOF Then
            datProbe.Recordset.Edit
            datProbe.Recordset.Update
        End If
    End If
    
End Sub



Private Sub fileImages_Click()

    If sActiveGrid = "Stimuli" Then
        gridStimuli.Columns(1).Text = fileImages.FileName
    ElseIf sActiveGrid = "Probe" Then
        gridProbe.Columns(0).Text = fileImages.FileName
    End If
    
End Sub

Private Sub Form_Activate()
    
    If gbRecalculateJuicePellets = True Then
        
        bLoadingRecords = True
        bLoadingForm = True
        
        gridStimuli.Visible = False
        gridProbe.Visible = False
        gridTrialInfo.Visible = False
        frameStatus.Caption = ""
        frameStatus.ZOrder
        frameStatus.Visible = True
        
        datSessionInfo.Refresh
        
        subClearApparatus
        subAssignApparatus
        subCountJuice
        subCountPellets
        subCalcApparatusPercentages
        
        frameStatus.Visible = False
        gridStimuli.Visible = True
        gridProbe.Visible = True
        gridTrialInfo.Visible = True
        
        bLoadingForm = False
        bLoadingRecords = False
    
    Else
    
        gsubScreenDimensions Me, 0, 0
        
        bLoadingForm = True
        bLoadingRecords = True
        bPrinting = False
        mnuClose.Enabled = False
        mnuCancelSetup.Enabled = True
        
        Refresh
        
        datSessionInfo.DatabaseName = gsProfile
        datSessionInfo.RecordSource = "SessionInfo"
        datSessionInfo.Refresh
    
        datTrialInfo.DatabaseName = gsProfile
        datTrialInfo.RecordSource = "TrialInfo"
        datTrialInfo.Refresh
    
        datStimuli.DatabaseName = gsProfile
        datStimuli.RecordSource = "StimuliInfo"
        datStimuli.Refresh
    
        datProbe.DatabaseName = gsProfile
        datProbe.RecordSource = "ProbeInfo"
        datProbe.Refresh
        
        datView.DatabaseName = gsProfile
        datView.RecordSource = "ViewConstants"
        datView.Refresh
        
        Refresh
        
        fileImages.Path = datSessionInfo.Recordset("ImagePath")
        
        If datSessionInfo.Recordset("DefaultVariableListLength") = True Then
            intListLength.Enabled = True
            intListLength.FillColor = COLOR_WHITE
        Else
            intListLength.Enabled = False
            intListLength.FillColor = COLOR_GRAY
        End If
        
        subLoadTrialRecords
        
        'may be changed during subLoadTrialRecords
        If bLoadingRecords = True Then
    
            subCountSames
            subCountDiffs
            subCountReinforces
            subCountJuice
            subCountPellets
            
            subCalcTypePercentages
            subCalcReinforcePercentages
            subCalcApparatusPercentages
            subCalculateStatistics
        
            datStimuli.RecordSource = "SELECT * FROM StimuliInfo WHERE TrialNum = 1"
            datStimuli.Refresh
            datProbe.RecordSource = "SELECT * FROM ProbeInfo WHERE TrialNum = 1"
            datProbe.Refresh
        
        End If
        
        bLoadingForm = False
        bLoadingRecords = False
        mnuClose.Enabled = True
        mnuCancelSetup.Enabled = False
    
    End If
    
End Sub

Private Sub Form_Load()

    gbRecalculateJuicePellets = False

    bLoadingForm = True
    
    sActiveGrid = "Stimuli"
   
End Sub

Private Sub gridProbe_Click()

    sActiveGrid = "Probe"
    
End Sub

Private Sub gridSames_HeadClick(ByVal ColIndex As Integer)

    subCalculateStatistics
    
End Sub


Private Sub gridSames_InitColumnProps()

    gridSames.Columns(0).Caption = "Config"
    gridSames.Columns(1).Caption = "Count"

End Sub

Private Sub gridStimuli_Click()

    sActiveGrid = "Stimuli"
    
End Sub

Private Sub gridTrialInfo_Change()

    Refresh
    
End Sub

Private Sub intCorrectPosition_Change()
    
    If bLoadingForm = False And bPrinting = False Then
    
        If CInt(intCorrectPosition.Text) > CInt(intListLength.Text) Then
            frameMessage.Visible = True
        Else
            frameMessage.Visible = False
        End If
    
    End If
        
End Sub

Private Sub intCorrectPosition_KeyUp(KeyCode As Integer, Shift As Integer)

    If bLoadingForm = False Then
        subCalculateStatistics
    End If

End Sub

Private Sub mnuCancelSetup_Click()
    
    Dim iRet As Integer
    
    If bLoadingRecords = True Then
        iRet = MsgBox("The trial data is not complete.  Do you to cancel the setup?", vbYesNo, "Setup Close")
        If iRet = vbYes Then
            bLoadingRecords = False
            mnuClose.Enabled = True
            mnuCancelSetup.Enabled = False
        End If
    End If

End Sub

Private Sub mnuClose_Click()
        
    Unload Me
    
End Sub

Private Sub subCalculateStatistics()
    
    datTrialInfo.Recordset.Edit
    datTrialInfo.Recordset("CorrectPosition") = intCorrectPosition
    datTrialInfo.Recordset.Update
    
    datStats.DatabaseName = gsProfile
    'datStats.RecordSource = "SELECT CorrectPosition, COUNT(*) AS Quantity FROM TrialInfo WHERE TrialType = 'Same' GROUP BY CorrectPosition"
    datStats.RecordSource = "SELECT Configuration, COUNT(*) AS Quantity FROM TrialInfo GROUP BY Configuration"
    datStats.Refresh

End Sub



Private Sub mnuPrint_Click()

    subPrintSetup

End Sub


Private Sub mnuUpdateProbabilities_Click()

    frmProbabilities.Show 1
    
End Sub
Public Sub subAssignStim()
    Dim iTrials As Integer
    Dim i As Integer
        
    Randomize
    iTrials = datSessionInfo.Recordset("NumOfTrials")
    datTrialInfo.Recordset.MoveFirst
    frameStatus.ZOrder
    frameStatus.Visible = True
    
    For i = 1 To iTrials
        datStimuli.RecordSource = "SELECT * FROM StimuliInfo WHERE TrialNum = " & datTrialInfo.Recordset("TrialNum")
        datStimuli.Refresh
        datStimuli.Recordset.Edit
        datProbe.RecordSource = "SELECT * FROM ProbeInfo WHERE TrialNum = " & datTrialInfo.Recordset("TrialNum")
        datProbe.Refresh
        datProbe.Recordset.Edit
        frameStatus.Caption = "Assigning Stimuli:  " & i & " of " & iTrials
        Refresh
        
        'size of stimuli array is >= display size + 1 (must account for probe item!)
        Dim colors(18) As Integer
        Dim set_size, display_size As Integer
        
        set_size = 18 'Int(datSessionInfo.Recordset("SetSize"))
        'the largest display we can have is
        '1 less than the size display size needs to equal 9 or less...
        display_size = 7   'set_size - 967
        
'        If (display_size > (set_size + 1)) Then
'            'MsgBox "Display Size >= set size. Program is going to produce unintended results.", vbCritical
'            display_size = set_size - 1
'        End If
        
        'grab the first 2 line segments randomly
        
        Dim c As Integer
        For c = 1 To display_size
            colors(c) = Int((set_size * Rnd) + 1)
        Next c
        
        'now we need to loop from 3 to display size and randomly pick new colors
        'a function that does sample without replacement is the best way to do this,
        'next best is just to have lots of nested loops
        Dim cc As Integer
        Dim go_on As Boolean
        'we already have the first stim (c=1) and we need to
        'account for the change-to color on the probe display
        c = display_size + 1
        Do
            colors(c) = Int((set_size * Rnd) + 1)
                
            'the logic here is that go_on = true only if NONE of the past colors
            'are equal to the current color. First occurence of go_on = false
            'breaks out of the loop.
            cc = 1
            Do
                go_on = colors(cc) <> colors(c)
                cc = cc + 1
            Loop While cc < c And go_on = True
        Loop While Not go_on
        
        'always going to be shape == 'c'
        Dim name As String
        Call intToColor(colors(1), name)
        datStimuli.Recordset("FileName") = name & ".bmp"
        
        Call intToColor(colors(2), name)
        datProbe.Recordset("Filename") = name & ".bmp"
        
        For c = 3 To display_size
            Call intToColor(colors(c), name)
            datProbe.Recordset("Filename" & (c - 1)) = name & ".bmp"
        Next c
        
        'fill the probe file name
        Call intToColor(colors(display_size + 1), name)
        datProbe.Recordset("ProbeFileName") = name & ".bmp"
       
        datStimuli.Recordset.Update
        datStimuli.Recordset.MoveNext
        datProbe.Recordset.Update
        datTrialInfo.Recordset.MoveNext
    Next i
        
    frameStatus.Visible = False
End Sub

Public Sub intToColor(ByVal Index As Integer, ByRef color_name As String)
    color_name = "s" & Index
End Sub

Public Sub PickWoReplaceD(Size As Integer, Pick)
      
    Dim Flag As Integer
    Dim CHECK As Integer
    
    Randomize
    Flag = 0
    For CHECK = 1 To Size
        If PickBuffD(CHECK) = 0 Then
            Flag = 1
        End If
    Next CHECK
    
    If Flag = 1 Then
AGAINd:
        Pick = Int(Rnd * Size) + 1
        If PickBuffD(Pick) = 1 Then
            GoTo AGAINd
        End If
        PickBuffD(Pick) = 1
    ElseIf Flag = 0 Then
        Erase PickBuffD
        Pick = Int(Rnd * Size) + 1
        PickBuffD(Pick) = 1
    End If

End Sub
Public Sub subAssignSerialPosition()
    Dim pos As Integer
    Dim iTrials As Integer
    Dim t As Integer
    Dim leng As Integer
    Dim Fi(10) As Integer
    Dim go As Boolean
    Dim amount As Integer
    Dim tem As Integer
    
    Randomize
    iTrials = datSessionInfo.Recordset("NumOfTrials")
    leng = datSessionInfo.Recordset("DefaultMaxListLength")
    amount = (iTrials / 2) / leng
    datTrialInfo.Recordset.MoveFirst
    For t = 1 To iTrials
        If datTrialInfo.Recordset("TrialType") = "Same" Then
            Do
                go = True
                pos = Int((leng * Rnd) + 1)
                Fi(pos) = Fi(pos) + 1
                If Fi(pos) > amount Then go = False
            Loop Until go = True
            datTrialInfo.Recordset.Edit
            datTrialInfo.Recordset("CorrectPosition") = pos
            datTrialInfo.Recordset.Update
        End If
        datTrialInfo.Recordset.MoveNext
    Next t

End Sub



Public Sub subAssignVariableProbeDelays()

    Dim pd As Integer 'probe delay
    Dim iTrials As Integer
    Dim t As Integer
    Dim DPD(8) As Integer 'DifferentProbeDelays
    Dim go As Boolean
    Dim amount As Integer
    
    
    datProbeInfo.DatabaseName = gsProfile
    datProbeInfo.RecordSource = "ProbeInfo"
    datProbeInfo.Refresh

    Randomize
    iTrials = datSessionInfo.Recordset("NumOfTrials")
    Dim randDelay As Boolean
    randDelay = False
    If (datSessionInfo.Recordset("VariableProbeDelay") = 0) Then
        randDelay = True
    End If
    
    amount = (iTrials) / 8
    datTrialInfo.Recordset.MoveFirst
    datProbeInfo.Recordset.MoveFirst
    For t = 1 To iTrials
         
        If (randDelay) Then
            Do
                go = True
                pd = Int((8 * Rnd) + 1)
                DPD(pd) = DPD(pd) + 1
                If DPD(pd) > amount Then go = False
            Loop Until go = True
        Else
            pd = Int((t - 1) / amount) + 1
        End If
        datProbeInfo.Recordset.Edit
        
        datProbeInfo.Recordset("Delay") = datView.Recordset("pd" & pd)
        
        datProbeInfo.Recordset.Update
        datProbeInfo.Recordset.MoveNext
        datTrialInfo.Recordset.MoveNext
    Next t

End Sub

