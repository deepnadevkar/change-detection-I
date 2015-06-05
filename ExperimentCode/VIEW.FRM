VERSION 5.00
Object = "{BDC217C8-ED16-11CD-956C-0000C04E4C0A}#1.1#0"; "TABCTL32.OCX"
Object = "{8D650141-6025-11D1-BC40-0000C042AEC0}#3.0#0"; "ssdw3b32.ocx"
Begin VB.Form frmViewData 
   BackColor       =   &H00808000&
   BorderStyle     =   1  'Fixed Single
   Caption         =   "View Data"
   ClientHeight    =   4140
   ClientLeft      =   1545
   ClientTop       =   1815
   ClientWidth     =   7665
   ControlBox      =   0   'False
   LinkTopic       =   "Form1"
   MaxButton       =   0   'False
   MinButton       =   0   'False
   PaletteMode     =   1  'UseZOrder
   ScaleHeight     =   4140
   ScaleWidth      =   7665
   Begin VB.Data datStimResponses 
      Caption         =   "StimResponses"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   3150
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   3720
      Visible         =   0   'False
      Width           =   2445
   End
   Begin VB.Data datTrialResults 
      Caption         =   "TrialResults"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   3150
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   3420
      Visible         =   0   'False
      Width           =   2445
   End
   Begin VB.Data datSessionResults 
      Caption         =   "SessionResults"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   3150
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   3120
      Visible         =   0   'False
      Width           =   2445
   End
   Begin VB.Data datProbeInfo 
      Caption         =   "ProbeInfo"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   3150
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   2820
      Visible         =   0   'False
      Width           =   2445
   End
   Begin VB.Data datStimuliInfo 
      Caption         =   "StimuliInfo"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   3150
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   2520
      Visible         =   0   'False
      Width           =   2445
   End
   Begin VB.Data datTrialInfo 
      Caption         =   "TrialInfo"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   3150
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   2220
      Visible         =   0   'False
      Width           =   2445
   End
   Begin VB.Data datSessionInfo 
      Caption         =   "SessionInfo"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   3150
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   1920
      Visible         =   0   'False
      Width           =   2445
   End
   Begin TabDlg.SSTab tabFileLevel 
      Height          =   6405
      Left            =   90
      TabIndex        =   0
      Top             =   120
      Width           =   9375
      _ExtentX        =   16536
      _ExtentY        =   11298
      _Version        =   393216
      Tabs            =   2
      TabsPerRow      =   2
      TabHeight       =   529
      BeginProperty Font {0BE35203-8F91-11CE-9DE3-00AA004BB851} 
         Name            =   "MS Sans Serif"
         Size            =   8.25
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      TabCaption(0)   =   "Session File:"
      TabPicture(0)   =   "VIEW.frx":0000
      Tab(0).ControlEnabled=   -1  'True
      Tab(0).Control(0)=   "tabSession"
      Tab(0).Control(0).Enabled=   0   'False
      Tab(0).ControlCount=   1
      TabCaption(1)   =   "Results File:"
      TabPicture(1)   =   "VIEW.frx":001C
      Tab(1).ControlEnabled=   0   'False
      Tab(1).Control(0)=   "tabResults"
      Tab(1).ControlCount=   1
      Begin TabDlg.SSTab tabSession 
         Height          =   5415
         Left            =   270
         TabIndex        =   1
         Top             =   750
         Width           =   8835
         _ExtentX        =   15584
         _ExtentY        =   9551
         _Version        =   393216
         Tabs            =   4
         Tab             =   2
         TabsPerRow      =   4
         TabHeight       =   529
         BackColor       =   12632256
         BeginProperty Font {0BE35203-8F91-11CE-9DE3-00AA004BB851} 
            Name            =   "MS Sans Serif"
            Size            =   8.25
            Charset         =   0
            Weight          =   400
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         TabCaption(0)   =   "Session Info"
         TabPicture(0)   =   "VIEW.frx":0038
         Tab(0).ControlEnabled=   0   'False
         Tab(0).Control(0)=   "gridSessionInfo"
         Tab(0).ControlCount=   1
         TabCaption(1)   =   "Trial Info"
         TabPicture(1)   =   "VIEW.frx":0054
         Tab(1).ControlEnabled=   0   'False
         Tab(1).Control(0)=   "gridTrialInfo"
         Tab(1).ControlCount=   1
         TabCaption(2)   =   "Stimuli Info"
         TabPicture(2)   =   "VIEW.frx":0070
         Tab(2).ControlEnabled=   -1  'True
         Tab(2).Control(0)=   "gridStimuliInfo"
         Tab(2).Control(0).Enabled=   0   'False
         Tab(2).ControlCount=   1
         TabCaption(3)   =   "Probe Info"
         TabPicture(3)   =   "VIEW.frx":008C
         Tab(3).ControlEnabled=   0   'False
         Tab(3).Control(0)=   "gridProbeInfo"
         Tab(3).ControlCount=   1
         Begin SSDataWidgets_B.SSDBGrid gridProbeInfo 
            Bindings        =   "VIEW.frx":00A8
            Height          =   4695
            Left            =   -74790
            TabIndex        =   8
            Top             =   540
            Width           =   8415
            _Version        =   196617
            ForeColorEven   =   0
            BackColorOdd    =   16777215
            RowHeight       =   423
            Columns(0).Width=   3200
            UseDefaults     =   0   'False
            _ExtentX        =   14843
            _ExtentY        =   8281
            _StockProps     =   79
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
         Begin SSDataWidgets_B.SSDBGrid gridStimuliInfo 
            Bindings        =   "VIEW.frx":00C3
            Height          =   4695
            Left            =   210
            TabIndex        =   5
            Top             =   540
            Width           =   8415
            _Version        =   196617
            ForeColorEven   =   0
            BackColorOdd    =   16777215
            RowHeight       =   423
            Columns(0).Width=   3200
            UseDefaults     =   0   'False
            _ExtentX        =   14843
            _ExtentY        =   8281
            _StockProps     =   79
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
         Begin SSDataWidgets_B.SSDBGrid gridTrialInfo 
            Bindings        =   "VIEW.frx":00E0
            Height          =   4695
            Left            =   -74790
            TabIndex        =   4
            Top             =   540
            Width           =   8415
            _Version        =   196617
            ForeColorEven   =   0
            BackColorOdd    =   16777215
            RowHeight       =   423
            Columns(0).Width=   3200
            UseDefaults     =   0   'False
            _ExtentX        =   14843
            _ExtentY        =   8281
            _StockProps     =   79
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
         Begin SSDataWidgets_B.SSDBGrid gridSessionInfo 
            Bindings        =   "VIEW.frx":00FB
            Height          =   4695
            Left            =   -74790
            TabIndex        =   3
            Top             =   540
            Width           =   8415
            _Version        =   196617
            ForeColorEven   =   0
            BackColorOdd    =   16777215
            RowHeight       =   423
            Columns(0).Width=   3200
            UseDefaults     =   0   'False
            _ExtentX        =   14843
            _ExtentY        =   8281
            _StockProps     =   79
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
      End
      Begin TabDlg.SSTab tabResults 
         Height          =   5415
         Left            =   -74730
         TabIndex        =   2
         Top             =   750
         Width           =   8835
         _ExtentX        =   15584
         _ExtentY        =   9551
         _Version        =   393216
         TabHeight       =   529
         BackColor       =   12632256
         BeginProperty Font {0BE35203-8F91-11CE-9DE3-00AA004BB851} 
            Name            =   "MS Sans Serif"
            Size            =   8.25
            Charset         =   0
            Weight          =   400
            Underline       =   0   'False
            Italic          =   0   'False
            Strikethrough   =   0   'False
         EndProperty
         TabCaption(0)   =   "Session Results"
         Tab(0).ControlEnabled=   -1  'True
         Tab(0).Control(0)=   "gridSessionResults"
         Tab(0).Control(0).Enabled=   0   'False
         Tab(0).ControlCount=   1
         TabCaption(1)   =   "Trial Results"
         Tab(1).ControlEnabled=   0   'False
         Tab(1).Control(0)=   "gridTrialResults"
         Tab(1).ControlCount=   1
         TabCaption(2)   =   "Stim Responses"
         Tab(2).ControlEnabled=   0   'False
         Tab(2).Control(0)=   "gridStimResponses"
         Tab(2).ControlCount=   1
         Begin SSDataWidgets_B.SSDBGrid gridStimResponses 
            Bindings        =   "VIEW.frx":0118
            Height          =   4695
            Left            =   -74760
            TabIndex        =   9
            Top             =   540
            Width           =   8415
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
            ForeColorEven   =   0
            BackColorOdd    =   16777215
            RowHeight       =   423
            Columns(0).Width=   3200
            UseDefaults     =   0   'False
            _ExtentX        =   14843
            _ExtentY        =   8281
            _StockProps     =   79
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
            Bindings        =   "VIEW.frx":0137
            Height          =   4695
            Left            =   -74760
            TabIndex        =   7
            Top             =   540
            Width           =   8415
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
            ForeColorEven   =   0
            BackColorOdd    =   16777215
            RowHeight       =   423
            Columns(0).Width=   3200
            UseDefaults     =   0   'False
            _ExtentX        =   14843
            _ExtentY        =   8281
            _StockProps     =   79
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
         Begin SSDataWidgets_B.SSDBGrid gridSessionResults 
            Bindings        =   "VIEW.frx":0155
            Height          =   4695
            Left            =   240
            TabIndex        =   6
            Top             =   540
            Width           =   8415
            _Version        =   196617
            ForeColorEven   =   0
            BackColorOdd    =   16777215
            RowHeight       =   423
            Columns(0).Width=   3200
            UseDefaults     =   0   'False
            _ExtentX        =   14843
            _ExtentY        =   8281
            _StockProps     =   79
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
      End
   End
   Begin VB.Menu mnuFile 
      Caption         =   "&File"
      Begin VB.Menu mnuClose 
         Caption         =   "&Close"
      End
   End
End
Attribute VB_Name = "frmViewData"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False
Option Explicit

Private Sub Form_Activate()

    Refresh
    
    If gsProfile <> "" Then

        datSessionInfo.DatabaseName = gsProfile
        datSessionInfo.RecordSource = "SessionInfo"
        datSessionInfo.Refresh
    
        datTrialInfo.DatabaseName = gsProfile
        datTrialInfo.RecordSource = "TrialInfo"
        datTrialInfo.Refresh
    
        datStimuliInfo.DatabaseName = gsProfile
        datStimuliInfo.RecordSource = "StimuliInfo"
        datStimuliInfo.Refresh
        
        datProbeInfo.DatabaseName = gsProfile
        datProbeInfo.RecordSource = "ProbeInfo"
        datProbeInfo.Refresh
    
    Else
    
        tabFileLevel.TabVisible(0) = False
    
    End If
    
        
    If gsResultsFile <> "" Then
        
        datSessionResults.DatabaseName = gsResultsFile
        datSessionResults.RecordSource = "SessionResults"
        datSessionResults.Refresh
        
        datTrialResults.DatabaseName = gsResultsFile
        datTrialResults.RecordSource = "TrialResults"
        datTrialResults.Refresh
        
        datStimResponses.DatabaseName = gsResultsFile
        datStimResponses.RecordSource = "StimResponses"
        datStimResponses.Refresh
    
    Else
    
        tabFileLevel.TabVisible(1) = False
    
    End If
    
End Sub

Private Sub Form_Load()
    
    gsubScreenDimensions Me, 0, 0

End Sub


Private Sub mnuClose_Click()

    Unload Me
    
End Sub



