import sqlite3
import re
import os
import sys

def parse_varint(data, i):
    val = 0
    shift = 0
    while i < len(data):
        b = data[i]
        i += 1
        val |= (b & 0x7f) << shift
        if not (b & 0x80):
            break
        shift += 7
    return val, i

def get_user_text(pl):
    pos = pl.find(b'\x12')
    while pos != -1:
        length, next_i = parse_varint(pl, pos + 1)
        if 0 < length <= len(pl) - next_i:
            cand = pl[next_i:next_i+length]
            try:
                s = cand.decode('utf-8')
                s_strip = s.strip()
                if len(s_strip) > 0 and not re.match(r'^[0-9a-f\-]{36}$', s_strip):
                    if not s_strip.startswith('call_') and s_strip not in ['USER', 'MODEL', 'SYSTEM']:
                        return s
            except Exception:
                pass
        pos = pl.find(b'\x12', pos + 1)
    return None

def get_model_text(pl):
    candidates = []
    pos = pl.find(b'\x0a')
    while pos != -1:
        length, next_i = parse_varint(pl, pos + 1)
        if 20 < length <= len(pl) - next_i:
            cand = pl[next_i:next_i+length]
            try:
                s = cand.decode('utf-8')
                s_strip = s.strip()
                if not s_strip.startswith('{') and not s_strip.startswith('call_') and 'sessionID' not in s_strip:
                    if len(s_strip) > 50:
                        candidates.append(s)
            except Exception:
                pass
        pos = pl.find(b'\x0a', pos + 1)
    if candidates:
        candidates.sort(key=len, reverse=True)
        return candidates[0]
    return None

def main():
    conv_id = "eb96079c-20ce-4f5e-bf42-ad8c2bb878f9"
    db_path = rf"C:\Users\Usuario\.gemini\antigravity\conversations\{conv_id}.db"
    output_md = "chat_export.md"
    
    if not os.path.exists(db_path):
        print(f"Database not found: {db_path}", file=sys.stderr)
        sys.exit(1)
        
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    
    # 1. Fetch user steps
    cur.execute("SELECT idx, step_payload FROM steps WHERE step_type = 14 ORDER BY idx;")
    user_steps = []
    for idx, pl in cur.fetchall():
        t = get_user_text(pl)
        if t:
            user_steps.append((idx, t))
            
    print(f"Total user turns detected: {len(user_steps)}")
    
    # 2. Match with model response
    turns = []
    for i, (u_idx, u_text) in enumerate(user_steps):
        next_u_idx = user_steps[i+1][0] if i + 1 < len(user_steps) else 9999999
        cur.execute(
            "SELECT idx, step_payload FROM steps WHERE step_type = 15 AND idx > ? AND idx < ? ORDER BY idx DESC;",
            (u_idx, next_u_idx)
        )
        rows = cur.fetchall()
        model_text = None
        m_idx = None
        for r_idx, pl in rows:
            mt = get_model_text(pl)
            if mt:
                model_text = mt
                m_idx = r_idx
                break
        turns.append({
            "turn_number": i + 1,
            "user_idx": u_idx,
            "user_text": u_text,
            "model_idx": m_idx,
            "model_text": model_text
        })
        
    # Write to Markdown
    with open(output_md, "w", encoding="utf-8") as f:
        f.write(f"# Transcripción Completa del Chat - Google Antigravity\n\n")
        f.write(f"- **ID de Conversación:** `{conv_id}`\n")
        f.write(f"- **Total de Turnos:** {len(turns)}\n")
        f.write(f"- **Generado automáticamente:** Sí\n\n")
        f.write(f"---\n\n")
        
        for turn in turns:
            num = turn["turn_number"]
            u_text = turn["user_text"].strip()
            # Clean possible XML tags
            u_text = re.sub(r'</?(?:USER_REQUEST|ADDITIONAL_METADATA)>', '', u_text).strip()
            
            f.write(f"## Turno {num}\n\n")
            f.write(f"### 👤 Usuario (Paso {turn['user_idx']}):\n\n")
            f.write(f"{u_text}\n\n")
            
            if turn["model_text"]:
                m_text = turn["model_text"].strip()
                f.write(f"### 🤖 Asistente (Antigravity - Paso {turn['model_idx']}):\n\n")
                f.write(f"{m_text}\n\n")
            else:
                f.write(f"### 🤖 Asistente (Antigravity):\n\n*(Turno en ejecución o sin respuesta de texto final)*\n\n")
                
            f.write(f"---\n\n")
            
    print(f"Exportación exitosa generada en: {os.path.abspath(output_md)}")
    print(f"Tamaño del archivo: {os.path.getsize(output_md)} bytes")

if __name__ == "__main__":
    main()

