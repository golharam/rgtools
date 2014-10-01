package tools.vcfbrowser;

import java.io.*;
import java.util.*;

public class AppSettings implements Serializable {

    static final long serialVersionUID = 4347888265825846218L;

    public AppSettings() {
    }

    public void setMRU(int num, String entry) {
        m_mruList[num] = entry;
    }

    public String getMRU(int num) {
        return m_mruList[num];
    }

    public void setMRUList(String[] mruList) {
        m_mruList = mruList;
    }

    public String[] getMRUList() {
        return m_mruList;
    }

    public void setUseExternalPhred(boolean b) {
        m_useExternalPhred = b;
    }

    public boolean getUseExternalPhred() {
        return m_useExternalPhred;
    }

    // deprecated
    public void setVectorSequence(String vecseq) {
        m_VectorSequence = vecseq;
    }

    // deprecated
    public String getVectorSequence() {
        return m_VectorSequence;
    }

    public void setVectorDatabase(String vecdb) {
        m_VectorDatabase = vecdb;
    }
    
    public String getVectorDatabase() {
        return m_VectorDatabase;
    }

    public LinkedList getVectorList() {
        if (m_Vectors == null)
            m_Vectors = new LinkedList();
            
        return m_Vectors;
    }
    
    public void setVectorList(LinkedList vectors) {
        m_Vectors = vectors;
    }
    
    public float getPhredCutoff() {
        return m_phredCutoff;
    }
    
    public void setPhredCutoff(float cutoff) {
        m_phredCutoff = cutoff;
    }
    
    private String[] m_mruList;
    private boolean m_useExternalPhred = false;
    private String m_VectorSequence; // deprecated
    private String m_VectorDatabase;
    private LinkedList m_Vectors;
    private float m_phredCutoff = (float)0.05;
}