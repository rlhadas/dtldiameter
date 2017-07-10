#include <windows.h>
#include "tree.cpp"

const char g_szClassName[] = "myWindowClass";

DTLInstance dtl;
std::string fileName = "TreeLifeData\\atest3.newick";//"TreeLifeData\\COG0500.newick"

int m_maxVScroll = 1000, m_maxHScroll = 1000;
int m_curVScroll = 0, m_curHScroll = 0;
const int m_w = 35, m_h = 50;

void loadTreeInfo() {
    Parser parser{};
    dtl = parser.parseFile(fileName);
	dtl.printInfo();
    string csvOut;
    dtl.computeAllStatistics(1,1,0,csvOut,false);
	dtl.printReconciliationGraph();
    cout << csvOut << endl;

    m_maxHScroll = 1200 + m_w*dtl.hostTree->size + m_w*dtl.paraTree->size 
                    + (dtl.hostTree->size+1)*(dtl.paraTree->size+1)*m_w/2;
    m_maxVScroll = 1200 + max<int>((dtl.hostTree->height+1)*(dtl.paraTree->height+1)*m_h,
                        (dtl.hostTree->height+1)*m_h*2*6);
}

void drawTree(HDC hdc, Tree* t, int x, int y, int w, int h) {
    if (t == nullptr || t->size < 1 || t->height < 1) return; // bad tree
    for (Vertex* v: t->preOrderVertexList) { 
        v->posX = x+(v->inOrderIndex*w)/(t->size);
        v->posY = y+h-(v->height*h)/(t->height);
    }
    for (Vertex* v: t->preOrderVertexList) { 
        int bw = 16, bh = 7; // half of the height/width of boxes
        RECT rect; rect.left  = v->posX-bw; rect.top    = v->posY-bh;
                   rect.right = v->posX+bw; rect.bottom = v->posY+bh;
        Rectangle(hdc, rect.left, rect.top, rect.right, rect.bottom);
        DrawText(hdc, v->label.c_str(), -1, &rect, DT_SINGLELINE|DT_CENTER);
        //TextOut(hdc, v->posX, v->posY, v->label.c_str(), v->label.length());
        if (v->left != nullptr) {
            MoveToEx(hdc, v->left->posX,  v->left->posY-bh, NULL);
            LineTo(  hdc, v->posX,        v->posY+bh);
            LineTo(  hdc, v->right->posX, v->right->posY-bh);
        }
    }
}

void drawReconciliationGraph(HDC hdc, vector<vector<MappingNode*>> mappingTable,
                                Tree* paraTree, Tree* hostTree,
                                int x, int y, int w, int h) {
    if (paraTree == nullptr || hostTree == nullptr || mappingTable.size() < 1) 
        return; // Not ready yet

    int groupW = (2*w)/(paraTree->size+1), groupH = h/(paraTree->height+1);
    for (Vertex* v: paraTree->preOrderVertexList) { 
        v->posX = x+(v->inOrderIndex*w)/(paraTree->size);
        v->posY = y+h-((v->height+1)*h)/(paraTree->height+1);
        Rectangle(hdc, v->posX, v->posY, v->posX+groupW, v->posY+groupH);
    }
    for (Vertex* v: hostTree->preOrderVertexList) { 
        v->posX = ((1+2*v->inOrderIndex)*groupW)/(2*hostTree->size);
        v->posY = groupH-((v->height+1)*groupH)/(hostTree->height+1);
    }
    for (Vertex* vp: paraTree->vertexList) { 
        for (Vertex* vh: hostTree->vertexList) { 
            MappingNode* m = mappingTable[vp->index][vh->index];
            if (!m->hit) continue;
            int tx = vp->posX+vh->posX, ty = vp->posY+vh->posY;
            m->posX = tx; m->posY = ty;
            int ew = 16, eh = 8; // half of the height/width of ellipse
            RECT rect; rect.left  = tx-ew; rect.top    = ty-eh;
                       rect.right = tx+ew; rect.bottom = ty+eh;
            Ellipse(hdc, rect.left, rect.top, rect.right, rect.bottom);
            DrawText(hdc, m->cstr().c_str(), -1, &rect, DT_SINGLELINE|DT_CENTER);
            //TextOut(hdc, vp->posX+vh->posX, vp->posY+vh->posY, 
            //        m->cstr().c_str(), m->cstr().length());
            int eventCnt = 0;
            for (EventNode* e: m->children) {
                int bw = 22, bh = 8; // half of the height/width of box
                e->posX = tx-(bw-ew)/2;// + (1+2*eventCnt++)*bw - (m->children.size())*bw;
                e->posY = ty + bh*2 + 2*eventCnt++*bh;
                RECT rec2; rec2.left  = e->posX-bw; rec2.top    = e->posY-bh;
                           rec2.right = e->posX+bw; rec2.bottom = e->posY+bh;
                Rectangle(hdc, rec2.left, rec2.top, rec2.right, rec2.bottom);
                DrawText(hdc, e->cstr().c_str(), -1, &rec2, DT_SINGLELINE|DT_CENTER);

                if (e->right != nullptr) {
                    MoveToEx(hdc, e->posX+bw,  e->posY, NULL);
                    LineTo(  hdc, e->right->posX, e->right->posY-eh);
                }
                if (e->left != nullptr) {
                    MoveToEx(hdc, e->posX-bw,  e->posY, NULL);
                    LineTo(  hdc, e->left->posX, e->left->posY-eh);
                }
            }
        }
    }

}

void drawReconciliationTree(HDC hdc, Tree* t, Tree* paraTree, ReconciliationGraph mappingTable,
                    ReconciliationTree rt, int x, int y, int w, int h, HPEN pen2) {
    if (t == nullptr || t->size < 1 || t->height < 1) return; // bad tree
    for (Vertex* v: t->preOrderVertexList) { 
        v->posX = x+(v->inOrderIndex*w)/(t->size);
        v->posY = y+h-(v->height*h)/(t->height);
    }
    int branchW = w/t->size, branchH = h/t->height;
    //cout << SSTR(branchW) << " " << SSTR(branchH) << endl;
    for (Vertex* v: paraTree->preOrderVertexList) { 
        v->posX = (((int)v->inOrderIndex - (int)paraTree->size/2)*branchW)/((int)paraTree->size);
        v->posY = -branchH/5-(v->height)*3*branchH/(5*paraTree->height);
    }
    for (Vertex* v: t->preOrderVertexList) { 
        int bw = 16, bh = 7; // half of the height/width of boxes
        RECT rect; rect.left  = v->posX-bw; rect.top    = v->posY-bh;
                   rect.right = v->posX+bw; rect.bottom = v->posY+bh;
        //Rectangle(hdc, rect.left, rect.top, rect.right, rect.bottom);
        DrawText(hdc, v->label.c_str(), -1, &rect, DT_SINGLELINE|DT_CENTER);
        //TextOut(hdc, v->posX, v->posY, v->label.c_str(), v->label.length());
        if (v->left != nullptr) {
            MoveToEx(hdc, v->left->posX,  v->left->posY-bh, NULL);
            LineTo(  hdc, v->posX,        v->posY+bh);
            LineTo(  hdc, v->right->posX, v->right->posY-bh);
        }
        if (v->isRoot()) {
            MoveToEx(hdc, v->posX,  v->posY-bh, NULL);
            LineTo(  hdc, v->posX,  v->posY-branchH);
        }
    }
    SelectObject(hdc, pen2);
    for (Vertex* vp: paraTree->vertexList) { 
        for (Vertex* vh: t->vertexList) { 
            MappingNode* m = mappingTable[vp->index][vh->index];
            if (!m->hit) continue;
            m->posX = vh->posX + vp->posX;
            m->posY = vh->posY + vp->posY;
            if (!vh->isRoot()) {
                m->posX += (vh->parent->posX - vh->posX)/(t->height + 7) * (vp->height + 2);
            }
            //TextOut(hdc, m->posX, m->posY, m->cstr().c_str(), m->cstr().length());
        }
    }
    for (EventNode* e: rt) {
        MappingNode* mp = e->parent;
        TextOut(hdc, mp->posX, mp->posY, e->cstr().c_str(), e->cstr().length());
        if (e->left != nullptr) {
            MoveToEx(hdc, mp->posX,  mp->posY, NULL);
            LineTo(hdc, e->left->posX,  e->left->posY);
        }
        if (e->right != nullptr) {
            MoveToEx(hdc, mp->posX,  mp->posY, NULL);
            LineTo(hdc, e->right->posX,  e->right->posY);
        }
    }
}



void drawDTL(HWND& hwnd) {
    PAINTSTRUCT ps;
    HDC hdc = BeginPaint (hwnd, &ps);
	SelectObject(hdc, GetStockObject(DEFAULT_GUI_FONT));
	SetBkMode(hdc, TRANSPARENT);

    HPEN hBluePen  = CreatePen(PS_SOLID, 1, RGB(0, 0, 255));
    HPEN hRedPen   = CreatePen(PS_SOLID, 1, RGB(200, 0, 0));
    HPEN hGreenPen = CreatePen(PS_SOLID, 1, RGB(0, 150, 0));
	//SetTextColor(hdc, RGB(255, 0, 0));

    int h = -m_curHScroll, v = -m_curVScroll;
    int th = 20, tv = 20;

    if (!(dtl.paraTree == nullptr || dtl.hostTree == nullptr || dtl.mappingTable.size() < 1)) { 
            
    SelectObject(hdc, hBluePen);
    TextOut(hdc, h+th-15, v+tv-20, "Host", 4);
    drawTree(hdc, dtl.hostTree, h+th, v+tv, m_w*dtl.hostTree->size, m_h*dtl.hostTree->height);
    th += m_w*dtl.hostTree->size + 40; 

    SelectObject(hdc, hRedPen);
    TextOut(hdc, h+th-15, v+tv-20, "Parasite", 8);
    drawTree(hdc, dtl.paraTree, h+th, v+tv, m_w*dtl.paraTree->size, m_h*dtl.paraTree->height);
    th += m_w*dtl.paraTree->size + 40; 

    SelectObject(hdc, hGreenPen);
    TextOut(hdc, h+th-15, v+tv-20, "RecGraph", 8);
    drawReconciliationGraph(hdc, dtl.mappingTable, dtl.paraTree, dtl.hostTree, 
                            h+th, v+tv, (dtl.hostTree->size+1)*(dtl.paraTree->size+1)*m_w/2, (dtl.hostTree->height+1)*(dtl.paraTree->height+1)*m_h);

    tv += 150 + max<int>(m_h*dtl.hostTree->height, m_h*dtl.paraTree->height);
    th = 20;

    SelectObject(hdc, hBluePen);
    TextOut(hdc, h+th-15, v+tv-20, "asymmetric median", 17);
    drawReconciliationTree(hdc, dtl.hostTree, dtl.paraTree, dtl.mappingTable,
                           dtl.rtAM, h+th, v+tv, 2*m_w*dtl.hostTree->size, 2*m_h*dtl.hostTree->height, hRedPen);
    tv += 150 + 2*(m_h*dtl.hostTree->height);

    SelectObject(hdc, hBluePen);
    TextOut(hdc, h+th-15, v+tv-20, "symmetric median", 16);
    drawReconciliationTree(hdc, dtl.hostTree, dtl.paraTree, dtl.mappingTable,
                           dtl.rtSM, h+th, v+tv, 2*m_w*dtl.hostTree->size, 2*m_h*dtl.hostTree->height, hRedPen);
    tv += 150 + 2*(m_h*dtl.hostTree->height);

    }
            
    DeleteObject(hBluePen);
    DeleteObject(hRedPen);
    DeleteObject(hGreenPen);
    EndPaint (hwnd, &ps);
}


void m_WM_VSCROLL(HWND hWnd, WPARAM wParam, LPARAM /*lParam*/)
{
    int nScrollCode = (int)LOWORD(wParam);
    int nPos = (short int)HIWORD(wParam);
    SCROLLINFO si = {sizeof(SCROLLINFO), 
                     SIF_PAGE|SIF_POS|SIF_RANGE|SIF_TRACKPOS, 0, 0, 0, 0, 
                     0};
    GetScrollInfo (hWnd, SB_VERT, &si);

    switch (nScrollCode)
    {
        case SB_THUMBPOSITION:
        case SB_THUMBTRACK:
          m_curVScroll = nPos;
          break;
        case SB_LINEDOWN:
            m_curVScroll += 100;
        break;
        case SB_LINEUP:
            m_curVScroll -= 100;
        break;
    }

    RECT r;
    GetWindowRect(hWnd, &r);

    si.fMask = SIF_POS | SIF_RANGE | SIF_PAGE;
    si.nMin = 0;
    si.nMax = m_maxVScroll;
    si.nPos = m_curVScroll;
    si.nPage = min<int>(si.nMax-1, r.bottom - r.top);

    SetScrollInfo (hWnd, SB_VERT, &si, TRUE);
}

void m_WM_HSCROLL(HWND hWnd, WPARAM wParam, LPARAM /*lParam*/)
{
    int nScrollCode = (int)LOWORD(wParam);
    int nPos = (short int)HIWORD(wParam);
    SCROLLINFO si = {sizeof(SCROLLINFO), 
                     SIF_PAGE|SIF_POS|SIF_RANGE|SIF_TRACKPOS, 0, 0, 0, 0, 
                     0};
    GetScrollInfo (hWnd, SB_HORZ, &si);

    switch (nScrollCode)
    {
        case SB_THUMBPOSITION:
        case SB_THUMBTRACK:
          m_curHScroll = nPos;
          break;
        case SB_LINEDOWN:
            m_curHScroll += 100;
        break;
        case SB_LINEUP:
            m_curHScroll -= 100;
        break;
    }

    RECT r;
    GetWindowRect(hWnd, &r);

    si.fMask = SIF_POS | SIF_RANGE | SIF_PAGE;
    si.nMin = 0;
    si.nMax = m_maxHScroll;
    si.nPos = m_curHScroll;
    si.nPage = min<int>(si.nMax-1, r.right - r.left);

    SetScrollInfo (hWnd, SB_HORZ, &si, TRUE);
}

// Step 4: the Window Procedure
LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
    switch(msg)
    {
        case WM_PAINT: 
            drawDTL(hwnd);
        break;
        case WM_CLOSE:
            DestroyWindow(hwnd);
        break;
        case WM_DESTROY:
            PostQuitMessage(0);
        break;
        case WM_VSCROLL:
            m_WM_VSCROLL(hwnd, wParam, lParam);
            RedrawWindow(hwnd, NULL, NULL, RDW_INVALIDATE | RDW_UPDATENOW | RDW_ERASE);
        break;
        case WM_HSCROLL:
            m_WM_HSCROLL(hwnd, wParam, lParam);
            RedrawWindow(hwnd, NULL, NULL, RDW_INVALIDATE | RDW_UPDATENOW | RDW_ERASE);
        break;
        default:
            return DefWindowProc(hwnd, msg, wParam, lParam);
    }
    return 0;
}

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE /*hPrevInstance*/,
    LPSTR /*lpCmdLine*/, int nCmdShow)
{
    WNDCLASSEX wc;
    HWND hwnd;
    MSG Msg;

    //Step 1: Registering the Window Class
    wc.cbSize        = sizeof(WNDCLASSEX);
    wc.style         = 0;
    wc.lpfnWndProc   = WndProc;
    wc.cbClsExtra    = 0;
    wc.cbWndExtra    = 0;
    wc.hInstance     = hInstance;
    wc.hIcon         = LoadIcon(NULL, IDI_APPLICATION);
    wc.hCursor       = LoadCursor(NULL, IDC_ARROW);
    wc.hbrBackground = (HBRUSH)(COLOR_WINDOW+1);
    wc.lpszMenuName  = NULL;
    wc.lpszClassName = g_szClassName;
    wc.hIconSm       = LoadIcon(NULL, IDI_APPLICATION);

    if(!RegisterClassEx(&wc))
    {
        MessageBox(NULL, "Window Registration Failed!", "Error!",
            MB_ICONEXCLAMATION | MB_OK);
        return 0;
    }

    // Step 2: Creating the Window
    hwnd = CreateWindowEx(
        WS_EX_CLIENTEDGE,
        g_szClassName,
        "DTL Visualizer",
        WS_OVERLAPPEDWINDOW | WS_HSCROLL | WS_VSCROLL,
        CW_USEDEFAULT, CW_USEDEFAULT, 120, 25,
        NULL, NULL, hInstance, NULL);

    if(hwnd == NULL)
    {
        MessageBox(NULL, "Window Creation Failed!", "Error!",
            MB_ICONEXCLAMATION | MB_OK);
        return 0;
    }

    ShowWindow(hwnd, nCmdShow);
    UpdateWindow(hwnd);
    loadTreeInfo();

    // Step 3: The Message Loop
    while(GetMessage(&Msg, NULL, 0, 0) > 0)
    {
        TranslateMessage(&Msg);
        DispatchMessage(&Msg);
    }
    return Msg.wParam;
}