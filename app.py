import streamlit as st
from langchain_openai import OpenAIEmbeddings
from langchain_community.vectorstores import FAISS
from langchain_openai import ChatOpenAI
from langchain.chains import ConversationalRetrievalChain
from langchain.text_splitter import CharacterTextSplitter
from langchain_community.document_loaders import TextLoader, Docx2txtLoader
import os
from dotenv import load_dotenv
from pathlib import Path

# Load environment variables
load_dotenv()

# Add these constants at the top of the file
VECTOR_STORE_DIR = "vector_store"
if not os.path.exists(VECTOR_STORE_DIR):
    os.makedirs(VECTOR_STORE_DIR)

def initialize_session_state():
    """Initialize session state variables"""
    session_vars = [
        'messages',
        'vectorstore',
        'processed_files'  # Track which files have been processed
    ]
    for var in session_vars:
        if var not in st.session_state:
            st.session_state[var] = [] if var == 'messages' or var == 'processed_files' else None

def process_document(uploaded_file):
    """Process document and update vectorstore"""
    # Save uploaded file temporarily
    with open(uploaded_file.name, "wb") as f:
        f.write(uploaded_file.getbuffer())
    
    # Load document based on file type
    if uploaded_file.name.endswith('.txt'):
        loader = TextLoader(uploaded_file.name)
    elif uploaded_file.name.endswith('.docx'):
        loader = Docx2txtLoader(uploaded_file.name)
    else:
        raise ValueError("Unsupported file format")
    
    documents = loader.load()
    
    # Split documents into chunks
    text_splitter = CharacterTextSplitter(chunk_size=1000, chunk_overlap=200)
    texts = text_splitter.split_documents(documents)
    
    # Create or update vectorstore
    embeddings = OpenAIEmbeddings()
    
    if st.session_state.vectorstore is None:
        # Create new vectorstore if none exists
        vectorstore = FAISS.from_documents(texts, embeddings)
    else:
        # Add new texts to existing vectorstore
        st.session_state.vectorstore.add_documents(texts)
        vectorstore = st.session_state.vectorstore
    
    # Save vectorstore
    save_path = os.path.join(VECTOR_STORE_DIR, "combined_vectorstore.faiss")
    vectorstore.save_local(save_path)
    
    # Clean up temporary file
    os.remove(uploaded_file.name)
    
    # Add to processed files list
    if uploaded_file.name not in st.session_state.processed_files:
        st.session_state.processed_files.append(uploaded_file.name)
    
    return vectorstore

def load_vectorstore():
    """Load combined vectorstore from disk"""
    save_path = os.path.join(VECTOR_STORE_DIR, "combined_vectorstore.faiss")
    if os.path.exists(save_path):
        embeddings = OpenAIEmbeddings()
        return FAISS.load_local(
            save_path, 
            embeddings,
            allow_dangerous_deserialization=True  # Only use this if you trust the source of the vector store
        )
    return None

def get_conversation_chain(vectorstore):
    """Create conversation chain"""
    llm = ChatOpenAI(temperature=0.7)
    conversation_chain = ConversationalRetrievalChain.from_llm(
        llm=llm,
        retriever=vectorstore.as_retriever(),
        return_source_documents=True
    )
    return conversation_chain

def handle_user_input(user_question):
    """Handle user input and generate response"""
    if st.session_state.vectorstore:
        conversation = get_conversation_chain(st.session_state.vectorstore)
        response = conversation({
            "question": user_question,
            "chat_history": [(msg["role"], msg["content"]) for msg in st.session_state.messages]
        })
        
        st.session_state.messages.append({"role": "user", "content": user_question})
        st.session_state.messages.append({"role": "assistant", "content": response["answer"]})

def main():
    st.set_page_config(page_title="Chat with Documents", page_icon="📚")
    initialize_session_state()

    st.header("Chat with Your Documents 📚")
    
    # Load existing vectorstore if available
    if st.session_state.vectorstore is None:
        st.session_state.vectorstore = load_vectorstore()
    
    # Show processed files
    if st.session_state.processed_files:
        st.sidebar.header("Processed Documents")
        st.sidebar.write("Documents in knowledge base:")
        for file in st.session_state.processed_files:
            st.sidebar.text(f"• {file}")
    
    # File upload
    uploaded_file = st.file_uploader(
        "Upload your document (TXT or DOCX)",
        type=["txt", "docx"]
    )

    if uploaded_file:
        if uploaded_file.name in st.session_state.processed_files:
            st.warning("This file has already been processed!")
        else:
            with st.spinner("Processing document..."):
                st.session_state.vectorstore = process_document(uploaded_file)
            st.success("Document processed and added to knowledge base!")

    # Chat interface
    if st.session_state.vectorstore:
        user_question = st.text_input("Ask a question about your documents:")
        if user_question:
            handle_user_input(user_question)

        # Display chat messages
        for message in st.session_state.messages:
            with st.chat_message(message["role"]):
                st.write(message["content"])

if __name__ == "__main__":
    main() 